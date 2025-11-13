import pickle
import os
from typing import Dict, Any, Set

from rdflib import Namespace, URIRef
try:
    from rdflib.plugins.parsers.ttl import TurtleParser
except Exception:
    from rdflib.plugins.parsers.notation3 import TurtleParser
from rdflib.parser import create_input_source
from urllib.parse import urljoin


DEFAULT_TTL = "/home/pickleopear/AlphaGEM/data_available/MNXmnet.ttl"
OUT_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxr_left_right.pkl"
ONLY_SIDE_CSV = "/home/pickleopear/AlphaGEM/data_available/only_side.csv"
MISSING_STOICH_CSV = "/home/pickleopear/AlphaGEM/data_available/mnxr_missing_stoich.csv"


SCHEMA = Namespace("https://rdf.metanetx.org/schema/")


def _last_fragment(s: str) -> str:
    if "#" in s:
        return s.rsplit("#", 1)[-1]
    return s.rstrip("/").rsplit("/", 1)[-1]


def extract_all_mnxr_with_coef(file_path: str) -> Dict[str, Dict[str, Dict[str, list]]]:
    # 映射与缓存
    reaction_node_to_id: Dict[object, str] = {}
    reaction_left_nodes: Dict[object, Set[object]] = {}
    reaction_right_nodes: Dict[object, Set[object]] = {}
    side_node_to_chems: Dict[object, Set[str]] = {}
    side_node_to_comps: Dict[object, Set[str]] = {}
    side_node_to_coef: Dict[object, float] = {}
    # 记录 mnx:side 指到的节点集合（与 left/right 结构一致）
    reaction_side_nodes: Dict[object, Set[object]] = {}

    # 结果（每个 MNXR 可能有多个具体反应节点）
    accum: Dict[str, Dict[str, Dict[str, list]]] = {}

    def ensure_bucket(mnxr: str, rxn_node: object) -> Dict[str, list]:
        if mnxr not in accum:
            accum[mnxr] = {}
        rkey = str(rxn_node)
        if rkey not in accum[mnxr]:
            accum[mnxr][rkey] = {"left": [], "right": []}
        return accum[mnxr][rkey]

    def on_triple(triple):
        s, p, o = triple
        # 标记 MNXR 反应
        if p == SCHEMA.mnxr and isinstance(o, URIRef):
            mnxr_id = _last_fragment(str(o))
            if mnxr_id.startswith("MNXR"):
                reaction_node_to_id[s] = mnxr_id
            return
        # 左/右边
        if p == SCHEMA.left:
            reaction_left_nodes.setdefault(s, set()).add(o)
            return
        if p == SCHEMA.right:
            reaction_right_nodes.setdefault(s, set()).add(o)
            return
        # side（未分配左右的一侧集合）
        if getattr(SCHEMA, "side", None) is not None and p == SCHEMA.side:
            reaction_side_nodes.setdefault(s, set()).add(o)
            return
        # 侧节点属性：chem / comp / coef
        if p == SCHEMA.chem and isinstance(o, URIRef):
            side_node_to_chems.setdefault(s, set()).add(_last_fragment(str(o)))
            return
        if getattr(SCHEMA, "comp", None) is not None and p == SCHEMA.comp and isinstance(o, URIRef):
            side_node_to_comps.setdefault(s, set()).add(_last_fragment(str(o)))
            return
        # 系数：coef（数字文本）
        if getattr(SCHEMA, "coef", None) is not None and p == SCHEMA.coef:
            try:
                side_node_to_coef[s] = float(str(o))
            except Exception:
                pass
            return

    # 解析
    parser = TurtleParser()
    source = create_input_source(location=file_path)

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

    # 组装结果（将每个 side 节点的 chem/comp/coef 展开到对应反应的左右）
    only_side_rows = []
    missing_stoich_mnxr = []
    
    for rxn_node, mnxr in reaction_node_to_id.items():
        bucket = ensure_bucket(mnxr, rxn_node)
        for side_node in reaction_left_nodes.get(rxn_node, set()):
            chems = side_node_to_chems.get(side_node, set())
            comps = side_node_to_comps.get(side_node, set()) or {''}
            coef = side_node_to_coef.get(side_node, 1.0)
            for ch in chems:
                for cp in comps:
                    bucket["left"].append({"chem": ch, "comp": cp, "coef": float(coef)})
        for side_node in reaction_right_nodes.get(rxn_node, set()):
            chems = side_node_to_chems.get(side_node, set())
            comps = side_node_to_comps.get(side_node, set()) or {''}
            coef = side_node_to_coef.get(side_node, 1.0)
            for ch in chems:
                for cp in comps:
                    bucket["right"].append({"chem": ch, "comp": cp, "coef": float(coef)})

        # 仅有 side 而无 left/right 的情况，记录样例到 CSV
        if not bucket["left"] and not bucket["right"]:
            has_side = rxn_node in reaction_side_nodes
            for side_node in reaction_side_nodes.get(rxn_node, set()):
                chems = side_node_to_chems.get(side_node, set())
                comps = side_node_to_comps.get(side_node, set()) or {''}
                coef = side_node_to_coef.get(side_node, 1.0)
                for ch in chems:
                    for cp in comps:
                        only_side_rows.append({
                            "mnxr": mnxr,
                            "reaction_key": str(rxn_node),
                            "chem": ch,
                            "comp": cp,
                            "coef": coef,
                        })
            # 如果既无 left/right 也无 side，则记为完全缺失
            if not has_side:
                missing_stoich_mnxr.append(mnxr)

    return accum, only_side_rows, missing_stoich_mnxr


def main(ttl_path: str = DEFAULT_TTL, out_path: str = OUT_PKL) -> None:
    data, only_side_rows, missing_stoich_mnxr = extract_all_mnxr_with_coef(ttl_path)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved {len(data)} MNXR entries with coefficients to {out_path}")

    # 输出 only_side.csv 以便人工核查
    if only_side_rows:
        import csv
        with open(ONLY_SIDE_CSV, "w", newline="") as fcsv:
            writer = csv.DictWriter(fcsv, fieldnames=["mnxr", "reaction_key", "chem", "comp", "coef"])
            writer.writeheader()
            writer.writerows(only_side_rows)
        print(f"Also wrote {len(only_side_rows)} rows of 'only side' entries to {ONLY_SIDE_CSV}")
    
    # 输出完全缺失计量信息的 MNXR
    if missing_stoich_mnxr:
        import csv
        with open(MISSING_STOICH_CSV, "w", newline="") as fcsv:
            writer = csv.writer(fcsv)
            writer.writerow(["mnxr"])
            for mnxr in missing_stoich_mnxr:
                writer.writerow([mnxr])
        print(f"Also wrote {len(missing_stoich_mnxr)} MNXR reactions without any stoichiometry to {MISSING_STOICH_CSV}")


if __name__ == "__main__":
    main()


