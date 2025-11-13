from rdflib import Namespace, URIRef, RDFS
try:
    from rdflib.plugins.parsers.ttl import TurtleParser
except Exception:
    from rdflib.plugins.parsers.notation3 import TurtleParser
from rdflib.parser import create_input_source
from urllib.parse import urlparse
import pickle


SCHEMA = Namespace("https://rdf.metanetx.org/schema/")
RDFS_COMMENT = RDFS.comment
# 同时兼容 http 与 https 的 chem 命名空间
CHEM_NS_SET = {
    "https://rdf.metanetx.org/chem/",
    "http://rdf.metanetx.org/chem/",
    "https://identifiers.org/metanetx.chemical:",
}


class StopParsing(Exception):
    pass


def _last_fragment(uri_or_str) -> str:
    s = str(uri_or_str)
    if "#" in s:
        return s.rsplit("#", 1)[-1]
    return s.rstrip("/").rsplit("/", 1)[-1]


def _pred_local_name(p: URIRef) -> str:
    return _last_fragment(p).lower()


def _extract_db_id_from_uri(uri: str) -> tuple[str, str] | None:
    """
    尝试从跨库 URI 提取 (db, accession)。
    规则：
    - 若末段含冒号，如 CHEBI:12345，则 db=CHEBI, id=CHEBI:12345
    - 否则依据域名/路径猜测常见库
    """
    last = _last_fragment(uri)
    if ":" in last:
        db = last.split(":", 1)[0]
        return db.upper(), last
    host = urlparse(uri).netloc.lower()
    if not host:
        return None
    if "kegg" in host:
        return "KEGG", last
    if "chebi" in host:
        return "CHEBI", last
    if "pubchem" in host:
        return "PUBCHEM", last
    if "hmdb" in host:
        return "HMDB", last
    if "rhea" in host:
        return "RHEA", last
    if "metacyc" in host or "biocyc" in host:
        return "METACYC", last
    return None


def read_mnxref_chems(
    file_path: str,
    target_mnxm_ids: set[str] | None = None,
    max_records: int | None = None,
) -> dict[str, dict]:
    """
    流式解析 MNXref.ttl 中的 chem 条目，提取注释：name, formula, charge, smiles 及 xrefs。

    返回：{ MNXMxx: {"name": str, "formula": str, "charge": int|str, "smiles": str,
                    "xrefs": {db: [ids...]}, "raw_xrefs": [(pred, obj_str), ...]} }
    - 若提供 target_mnxm_ids，则仅提取其子集
    - 若提供 max_records，则在达到数量后提前停止解析
    """
    results: dict[str, dict] = {}
    pred_stats: dict[str, int] = {}

    def ensure_entry(mnxm_id: str) -> dict:
        if mnxm_id not in results:
            results[mnxm_id] = {
                "name": None,
                "description": None,
                "formula": None,
                "charge": None,
                "smiles": None,
                "xrefs": {},  # db -> set(ids)
                "raw_xrefs": [],
            }
        return results[mnxm_id]

    def add_xref(entry: dict, db: str, acc: str) -> None:
        x = entry["xrefs"].setdefault(db, set())
        x.add(acc)

    def on_triple(triple: tuple[object, object, object]) -> None:
        s, p, o = triple
        if not isinstance(s, URIRef):
            return
        su = str(s)
        if not any(su.startswith(ns) for ns in CHEM_NS_SET):
            return
        mnxm_id = _last_fragment(su)
        if target_mnxm_ids is not None and mnxm_id not in target_mnxm_ids:
            return

        entry = ensure_entry(mnxm_id)
        pl = _pred_local_name(p)
        pred_stats[pl] = pred_stats.get(pl, 0) + 1
        # 捕获核心属性（大小写/变体容错）
        # 优先用 comment 填 name
        if pl == "comment":
            if entry["name"] is None and hasattr(o, "toPython"):
                entry["name"] = str(o)
            return
        # 其次 label
        if pl in {"description", "comment"}:
            if entry["description"] is None and hasattr(o, "toPython"):
                entry["description"] = str(o)
            return
        if pl == "formula":
            if entry["formula"] is None and hasattr(o, "toPython"):
                entry["formula"] = str(o)
            return
        if pl == "charge":
            if entry["charge"] is None and hasattr(o, "toPython"):
                try:
                    entry["charge"] = int(str(o))
                except Exception:
                    entry["charge"] = str(o)
            return
        if pl in {"smiles", "smile"}:
            if entry["smiles"] is None and hasattr(o, "toPython"):
                entry["smiles"] = str(o)
            return

        # 跨库 ID：专用谓词或通用链接；支持对象为 URI 或文本
        XREF_PRED_SET = {
            # MNXref schema 针对 CHEM 的交叉引用谓词
            "chemxref", "chemrefer", "chemsource",
            # 通用 cross-ref 与同义链接
            "xref", "xrefs", "dbxref", "is", "seealso", "sameas",
            # 常见库专用谓词（若数据中直接使用）
            "chebi", "kegg", "hmdb", "pubchem", "rhea", "metacyc",
            "bigg", "seed", "inchi", "inchikey",
        }
        if pl in XREF_PRED_SET:
            # 记录原始 xref
            entry["raw_xrefs"].append((pl, str(o)))
            IGNORED_DBS = {"XREF", "CHEMXREF", "CHEMSOURCE", "CHEMREFER"}
            if isinstance(o, URIRef):
                dbid = _extract_db_id_from_uri(str(o))
                if dbid is not None:
                    if dbid[0] not in IGNORED_DBS:
                        add_xref(entry, dbid[0], dbid[1])
                else:
                    fallback_db = pl.upper()
                    if fallback_db not in IGNORED_DBS:
                        add_xref(entry, fallback_db, _last_fragment(o))
            else:
                val = str(o)
                if ":" in val:
                    db, acc = val.split(":", 1)
                    dbu = db.upper()
                    if dbu not in IGNORED_DBS:
                        add_xref(entry, dbu, f"{db}:{acc}")
                else:
                    # 无库名前缀，忽略
                    pass
            return

        # 提前停止：达到数量上限
        if max_records is not None and len(results) >= max_records:
            raise StopParsing()

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

    from urllib.parse import urljoin

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
    try:
        parser.parse(source, graph_like)
    except StopParsing:
        pass

    # xrefs set -> list（过滤空项）；不附加统计键，保持纯注释输出
    for mnxm, entry in results.items():
        entry["xrefs"] = {db: sorted(list(ids)) for db, ids in entry["xrefs"].items() if ids}

    # 新增：打印出现最多的谓词（只打印前20条）
    print("谓词出现频率排名前20:")
    for k, v in sorted(pred_stats.items(), key=lambda x: -x[1])[:20]:
        print(f"{k}: {v}")

    return results


if __name__ == "__main__":
    # 全量读取（可能耗时，流式执行）
    ttl_path = "../data_available/MNXref.ttl"
    ann = read_mnxref_chems(ttl_path)
    # 单独导出规范化的 xrefs
    xrefs_only = {k: v.get("xrefs", {}) for k, v in ann.items() if k.startswith("MNXM")}
    with open("../../data_available/mnxref_chem_annotations.pkl", "wb") as f:
        pickle.dump(ann, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open("../../data_available/mnxref_chem_xrefs_only.pkl", "wb") as f:
        pickle.dump(xrefs_only, f, protocol=pickle.HIGHEST_PROTOCOL)


