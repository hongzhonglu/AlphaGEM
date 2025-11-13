import os
import pickle
from typing import Dict, Any, Set

from rdflib import URIRef
try:
    from rdflib.plugins.parsers.ttl import TurtleParser
except Exception:
    from rdflib.plugins.parsers.notation3 import TurtleParser
from rdflib.parser import create_input_source
from urllib.parse import urljoin


REAC_XREFS_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxref_reac_xrefs.pkl"
RHEA_RDF = "/home/pickleopear/AlphaGEM/rhea/rhea.rdf"
OUT_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxr_reversibility.pkl"


def _last_fragment(uri_or_str) -> str:
    s = str(uri_or_str)
    if "#" in s:
        return s.rsplit("#", 1)[-1]
    return s.rstrip("/").rsplit("/", 1)[-1]


def _pred_local_name(p: URIRef) -> str:
    return _last_fragment(p).lower()


def _normalize_rhea_id(val: str) -> str | None:
    if not val:
        return None
    s = str(val)
    frag = _last_fragment(s)
    if ":" in frag:
        frag = frag.split(":", 1)[-1]
    # Rhea IDs通常是数字，如 12345，或 RHEA:12345
    return frag if frag.isdigit() else None


def read_mnxr_xrefs(path: str) -> Dict[str, Dict[str, list]]:
    if not os.path.exists(path):
        return {}
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data if isinstance(data, dict) else {}


def read_rhea_direction(rhea_rdf_path: str, target_rhea_ids: Set[str]) -> Dict[str, str]:
    """
    解析 Rhea RDF，提取目标 Rhea ID 的方向：both/L2R/R2L
    以关键词匹配方式兼容不同谓词命名（如 direction/equationDirection）。
    """
    if not os.path.exists(rhea_rdf_path) or not target_rhea_ids:
        return {}

    results: Dict[str, str] = {}

    def on_triple(triple):
        s, p, o = triple
        if not isinstance(s, URIRef):
            return
        rid = _normalize_rhea_id(str(s))
        if rid is None or rid not in target_rhea_ids:
            return
        pl = _pred_local_name(p)
        if "direction" not in pl:
            return
        val = str(o)
        low = val.lower()
        if any(k in low for k in ["bidirectional", "both", "reversible", "bidir"]):
            results[rid] = "both"
        elif "left" in low and "right" in low:
            # 如 left-to-right
            if "left-to-right" in low or "ltr" in low:
                results[rid] = "L2R"
            elif "right-to-left" in low or "rtl" in low:
                results[rid] = "R2L"
        elif "forward" in low:
            results[rid] = "L2R"
        elif "reverse" in low:
            results[rid] = "R2L"

    parser = TurtleParser()
    source = create_input_source(location=rhea_rdf_path)

    class _Store:
        context_aware = True

    class _NSMgr:
        def __init__(self):
            self._ns = {}
        def bind(self, prefix, uri):
            self._ns[prefix] = uri

    class _GraphLike:
        def __init__(self):
            self.store = _Store()
            self.namespace_manager = _NSMgr()
            self._base = source.getPublicId() or source.getSystemId() or ""
            self.base = self._base
        def bind(self, prefix, uri):
            self.namespace_manager.bind(prefix, uri)
        def add(self, triple):
            on_triple(triple)
        def addN(self, quads):
            for s, p, o, _ in quads:
                on_triple((s, p, o))
        def absolutize(self, uri, defrag=True):
            return urljoin(self._base, str(uri))
        def skolem_genid(self, bnode):
            return None

    graph_like = _GraphLike()
    parser.parse(source, graph_like)
    return results


def build_mnxr_reversibility_map() -> Dict[str, Dict[str, Any]]:
    mnxr_xrefs = read_mnxr_xrefs(REAC_XREFS_PKL)
    # 收集所有 Rhea ID
    rhea_ids: Set[str] = set()
    for mnxr, xmap in mnxr_xrefs.items():
        if not isinstance(xmap, dict):
            continue
        for db, ids in xmap.items():
            if db.upper() != "RHEA":
                continue
            for rid in ids or []:
                norm = _normalize_rhea_id(rid)
                if norm:
                    rhea_ids.add(norm)

    rhea_dir = read_rhea_direction(RHEA_RDF, rhea_ids)

    out: Dict[str, Dict[str, Any]] = {}
    for mnxr, xmap in mnxr_xrefs.items():
        dir_label = None
        if isinstance(xmap, dict):
            rids = xmap.get("RHEA") or xmap.get("rhea") or []
            for rid in rids:
                norm = _normalize_rhea_id(rid)
                if norm and norm in rhea_dir:
                    dir_label = rhea_dir[norm]
                    break
        if dir_label is None:
            reversible = None
        else:
            reversible = True if dir_label == "both" else False
        out[mnxr] = {"direction": dir_label, "reversible": reversible}
    return out


if __name__ == "__main__":
    revmap = build_mnxr_reversibility_map()
    with open(OUT_PKL, "wb") as f:
        pickle.dump(revmap, f)
    print(f"Saved reversibility info for {len(revmap)} MNXR reactions -> {OUT_PKL}")


