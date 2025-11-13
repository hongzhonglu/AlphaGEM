from rdflib import Namespace, URIRef, BNode 
try:
    from rdflib.plugins.parsers.ttl import TurtleParser
except Exception:
    from rdflib.plugins.parsers.notation3 import TurtleParser
from rdflib.parser import create_input_source
from urllib.parse import urlparse, urljoin
import pickle


SCHEMA = Namespace("https://rdf.metanetx.org/schema/")
# 反应主体命名空间（支持 https/http 以及 identifiers 代理）
REAC_SUBJ_PREFIXES = {
    "https://rdf.metanetx.org/reac/",
    "http://rdf.metanetx.org/reac/",
    "https://identifiers.org/metanetx.reaction:",
}
# 化学物命名空间（用于提取左右代谢物）
CHEM_NS_SET = {
    "https://rdf.metanetx.org/chem/",
    "http://rdf.metanetx.org/chem/",
    "https://identifiers.org/metanetx.chemical:",
}
# 区室命名空间
COMP_NS_SET = {
    "https://rdf.metanetx.org/comp/",
    "http://rdf.metanetx.org/comp/",
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
    last = _last_fragment(uri)
    if ":" in last:
        db = last.split(":", 1)[0]
        return db.upper(), last
    host = urlparse(uri).netloc.lower()
    if not host:
        return None
    if "kegg" in host:
        return "KEGG", last
    if "metacyc" in host or "biocyc" in host:
        return "METACYC", last
    if "rhea" in host:
        return "RHEA", last
    if "reactome" in host:
        return "REACTOME", last
    if "bigg" in host:
        return "BIGG", last
    if "seed" in host or "modelseed" in host:
        return "SEED", last
    if "expasy" in host or "enzyme" in host:
        return "EC", last
    return None


def read_mnxref_reacs(
    file_path: str,
    target_mnxr_ids: set[str] | None = None,
    max_records: int | None = None,
) -> tuple[dict[str, dict], list[str]]:
    """
    流式解析 MNXref.ttl 中的 reaction (MNXR*)，提取 xref 信息。

    返回：(annotations_dict, missing_mnxr_list)
    - annotations_dict: { MNXRxx: {"name": str|None, "xrefs": {db: [ids...]}, "raw_xrefs": [(pred, obj_str), ...], "ec": [..]} }
    - missing_mnxr_list: 在数据库中有但缺失详细信息的 MNXR ID
    """
    results: dict[str, dict] = {}
    all_mnxr_ids: set[str] = set()  # 收集所有出现过的 MNXR ID

    def ensure_entry(mnxr_id: str) -> dict:
        if mnxr_id not in results:
            results[mnxr_id] = {
                "name": None,
                "xrefs": {},      # db -> set(ids)
                "raw_xrefs": [],  # (pred, obj)
                "ec": set(),      # set of EC strings
            }
        return results[mnxr_id]

    def add_xref(entry: dict, db: str, acc: str) -> None:
        x = entry["xrefs"].setdefault(db, set())
        x.add(acc)

    def on_triple(triple: tuple[object, object, object]) -> None:
        s, p, o = triple
        if not isinstance(s, URIRef):
            return
        su = str(s)
        if not any(su.startswith(pref) for pref in REAC_SUBJ_PREFIXES):
            return
        mnxr_id = _last_fragment(su)
        if not mnxr_id.upper().startswith("MNXR"):
            return
        if target_mnxr_ids is not None and mnxr_id not in target_mnxr_ids:
            return

        # 记录所有出现的 MNXR ID
        all_mnxr_ids.add(mnxr_id)

        entry = ensure_entry(mnxr_id)
        pl = _pred_local_name(p)

        # EC 号（直接谓词或分类字段中可能出现）
        EC_PRED_SET = {"ec", "ecnumber", "ec_number", "ec-code"}
        if pl in EC_PRED_SET:
            if isinstance(o, URIRef):
                entry["ec"].add(_last_fragment(o).split(":", 1)[-1])
            else:
                entry["ec"].add(str(o).split(":", 1)[-1])
        if pl == "classification" and not isinstance(o, URIRef):
            txt = str(o).strip()
            # 简单匹配 a.b.c.d 形式
            import re
            if re.match(r"^\d+\.(\d+|n)\.(\d+|n)\.(\d+|n)$", txt, re.IGNORECASE):
                entry["ec"].add(txt.replace("EC ", "").replace("ec ", ""))

        # 名称/标签
        if pl in {"name", "label"}:
            if entry.get("name") is None and hasattr(o, "toPython"):
                entry["name"] = str(o)
        # 参考 MNXref schema: reacXref / reacRefer / reacSource
        XREF_PRED_SET = {
            "reacxref", "reacrefer", "reacsource",
            # 通用 cross-ref 与同义链接（数据中若出现）
            "xref", "xrefs", "dbxref", "is", "seealso", "sameas",
        }
        if pl in XREF_PRED_SET:
            entry["raw_xrefs"].append((pl, str(o)))
            IGNORED_DBS = {"XREF", "REACXREF", "REACSOURCE", "REACREFER"}
            if isinstance(o, URIRef):
                dbid = _extract_db_id_from_uri(str(o))
                if dbid is not None and dbid[0] not in IGNORED_DBS:
                    # 记录 xref
                    add_xref(entry, dbid[0], dbid[1])
                    # 若为 EC，抽取编号
                    if dbid[0] == "EC":
                        entry["ec"].add(dbid[1].split(":", 1)[-1])
                # 无法识别库名则忽略，避免引入占位类
            else:
                # 文本值：仅接受 DB:ACC 形式
                val = str(o)
                if ":" in val:
                    db, acc = val.split(":", 1)
                    dbu = db.upper()
                    if dbu not in IGNORED_DBS:
                        add_xref(entry, dbu, f"{db}:{acc}")
                        if dbu == "EC":
                            entry["ec"].add(acc)
            # 不 return，继续看是否有其它与反应相关的有用属性（目前仅收 xref）

        # 提前停止：达到数量上限（以反应实体计数）
        if max_records is not None and len(results) >= max_records:
            raise StopParsing()

    # 解析（流式）
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
    try:
        parser.parse(source, graph_like)
    except StopParsing:
        pass

    # xrefs set -> list；清理空项
    for mnxr, entry in results.items():
        entry["xrefs"] = {db: sorted(list(ids)) for db, ids in entry["xrefs"].items() if ids}
        entry["ec"] = sorted([e for e in entry.get("ec", set()) if e])
    
    # 识别缺失详细信息的 MNXR：有 ID 但无任何 name/xrefs/ec
    missing_mnxr = []
    for mnxr_id in all_mnxr_ids:
        entry = results.get(mnxr_id)
        if entry is None or (
            entry["name"] is None 
            and not entry["xrefs"] 
            and not entry["ec"] 
            and not entry["raw_xrefs"]
        ):
            missing_mnxr.append(mnxr_id)
    
    return results, sorted(missing_mnxr)


def read_mnxref_reacs_with_stoich(
    mnxref_path: str,
    mnxmnet_path: str,
    target_mnxr_ids: set[str] | None = None,
) -> tuple[dict[str, dict], dict[str, dict], list[str]]:
    """
    从 MNXref.ttl 提取注释信息，从 MNXref.ttl 和 MNXmnet.ttl 提取计量学信息
    （优先使用 MNXmnet.ttl 的计量学信息）
    
    返回：(annotations_dict, stoich_dict, missing_mnxr_list)
    - annotations_dict: { MNXRxx: {"name": ..., "xrefs": {...}, "ec": [...]} }
    - stoich_dict: { MNXRxx: {"reaction_nodes": {node: {"left": [...], "right": [...]}}}}
      其中 left/right 为 [{"chem": "MNXM...", "comp": "MNXC...", "coef": 1.0}, ...]
    - missing_mnxr_list: 缺失详细信息的 MNXR ID
    """
    # ===== 第一步：从 MNXref.ttl 提取注释信息和计量学信息 =====
    print("正在从 MNXref.ttl 提取注释信息和计量学信息...")
    annotations, stoich_from_ref, missing_ann = _extract_all_from_mnxref(mnxref_path, target_mnxr_ids)
    print(f"从 MNXref.ttl 提取了 {len(annotations)} 个反应的注释信息")
    print(f"从 MNXref.ttl 提取了 {len(stoich_from_ref)} 个反应的计量学信息")
    
    # ===== 第二步：从 MNXmnet.ttl 提取计量学信息（优先） =====
    print("正在从 MNXmnet.ttl 提取计量学信息...")
    stoich_from_mnet = _extract_stoich_from_mnxmnet(mnxmnet_path, target_mnxr_ids)
    print(f"从 MNXmnet.ttl 提取了 {len(stoich_from_mnet)} 个反应的计量学信息")
    
    # ===== 合并计量学信息：MNXmnet.ttl 优先，缺失的用 MNXref.ttl 补充 =====
    stoich_data = stoich_from_ref.copy()  # 先使用 MNXref.ttl 的数据
    for mnxr, mnet_data in stoich_from_mnet.items():
        if mnxr in stoich_data:
            # 如果 MNXmnet.ttl 有数据，优先使用（合并 reaction_nodes）
            for rxn_key, sides in mnet_data.get("reaction_nodes", {}).items():
                if rxn_key in stoich_data[mnxr]["reaction_nodes"]:
                    # 如果同一个 reaction_node，优先使用 MNXmnet.ttl 的数据
                    stoich_data[mnxr]["reaction_nodes"][rxn_key] = sides
                else:
                    stoich_data[mnxr]["reaction_nodes"][rxn_key] = sides
        else:
            # 如果 MNXref.ttl 没有，直接添加
            stoich_data[mnxr] = mnet_data
    
    print(f"合并后共有 {len(stoich_data)} 个反应的计量学信息")
    
    return annotations, stoich_data, missing_ann


def _extract_all_from_mnxref(
    file_path: str,
    target_mnxr_ids: set[str] | None = None,
) -> tuple[dict[str, dict], dict[str, dict], list[str]]:
    """
    从 MNXref.ttl 同时提取反应注释和计量学信息（通过空白节点）
    """
    # ===== 注释信息 =====
    annotations: dict[str, dict] = {}
    all_mnxr_ids: set[str] = set()
    
    def ensure_entry(mnxr_id: str) -> dict:
        if mnxr_id not in annotations:
            annotations[mnxr_id] = {
                "name": None,
                "xrefs": {},
                "raw_xrefs": [],
                "ec": set(),
            }
        return annotations[mnxr_id]
    
    def add_xref(entry: dict, db: str, acc: str) -> None:
        x = entry["xrefs"].setdefault(db, set())
        x.add(acc)
    
    # ===== 计量学信息 =====
    # 反应节点 -> MNXR ID
    reaction_node_to_id: dict[URIRef, str] = {}
    # 反应节点 -> left/right 空白节点
    reaction_left_bnodes: dict[URIRef, set[BNode]] = {}
    reaction_right_bnodes: dict[URIRef, set[BNode]] = {}
    # 空白节点 -> chem/comp/coef
    bnode_to_chems: dict[BNode, set[str]] = {}
    bnode_to_comps: dict[BNode, set[str]] = {}
    bnode_to_coef: dict[BNode, float] = {}
    
    def on_triple(triple: tuple[object, object, object]) -> None:
        s, p, o = triple
        pl = _pred_local_name(p)
        
        # ===== 处理反应节点的注释信息 =====
        if isinstance(s, URIRef):
            su = str(s)
            if any(su.startswith(pref) for pref in REAC_SUBJ_PREFIXES):
                mnxr_id = _last_fragment(su)
                if mnxr_id.upper().startswith("MNXR"):
                    if target_mnxr_ids is None or mnxr_id in target_mnxr_ids:
                        all_mnxr_ids.add(mnxr_id)
                        entry = ensure_entry(mnxr_id)
                        reaction_node_to_id[s] = mnxr_id
                        
                        # EC 号
                        EC_PRED_SET = {"ec", "ecnumber", "ec_number", "ec-code"}
                        if pl in EC_PRED_SET:
                            if isinstance(o, URIRef):
                                entry["ec"].add(_last_fragment(o).split(":", 1)[-1])
                            else:
                                entry["ec"].add(str(o).split(":", 1)[-1])
                        if pl == "classification" and not isinstance(o, URIRef):
                            txt = str(o).strip()
                            import re
                            if re.match(r"^\d+\.(\d+|n)\.(\d+|n)\.(\d+|n)$", txt, re.IGNORECASE):
                                entry["ec"].add(txt.replace("EC ", "").replace("ec ", ""))
                        
                        # 名称
                        if pl in {"name", "label"}:
                            if entry.get("name") is None and hasattr(o, "toPython"):
                                entry["name"] = str(o)
                        
                        # Xrefs
                        XREF_PRED_SET = {
                            "reacxref", "reacrefer", "reacsource",
                            "xref", "xrefs", "dbxref", "is", "seealso", "sameas",
                        }
                        if pl in XREF_PRED_SET:
                            entry["raw_xrefs"].append((pl, str(o)))
                            IGNORED_DBS = {"XREF", "REACXREF", "REACSOURCE", "REACREFER"}
                            if isinstance(o, URIRef):
                                dbid = _extract_db_id_from_uri(str(o))
                                if dbid is not None and dbid[0] not in IGNORED_DBS:
                                    add_xref(entry, dbid[0], dbid[1])
                                    if dbid[0] == "EC":
                                        entry["ec"].add(dbid[1].split(":", 1)[-1])
                            else:
                                val = str(o)
                                if ":" in val:
                                    db, acc = val.split(":", 1)
                                    dbu = db.upper()
                                    if dbu not in IGNORED_DBS:
                                        add_xref(entry, dbu, f"{db}:{acc}")
                                        if dbu == "EC":
                                            entry["ec"].add(acc)
                        
                        # ===== 处理反应的 left/right 空白节点 =====
                        if pl == "left" and isinstance(o, BNode):
                            reaction_left_bnodes.setdefault(s, set()).add(o)
                        elif pl == "right" and isinstance(o, BNode):
                            reaction_right_bnodes.setdefault(s, set()).add(o)
        
        # ===== 处理空白节点的 chem/comp/coef =====
        if isinstance(s, BNode):
            if pl == "chem" and isinstance(o, URIRef):
                ou = str(o)
                if any(ou.startswith(ns) for ns in CHEM_NS_SET):
                    bnode_to_chems.setdefault(s, set()).add(_last_fragment(ou))
            elif pl == "comp" and isinstance(o, URIRef):
                ou = str(o)
                if any(ou.startswith(ns) for ns in COMP_NS_SET):
                    bnode_to_comps.setdefault(s, set()).add(_last_fragment(ou))
            elif pl == "coef":
                try:
                    bnode_to_coef[s] = float(o)
                except Exception:
                    pass
    
    # 解析（流式）
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
    
    # ===== 清理注释信息 =====
    for mnxr, entry in annotations.items():
        entry["xrefs"] = {db: sorted(list(ids)) for db, ids in entry["xrefs"].items() if ids}
        entry["ec"] = sorted([e for e in entry.get("ec", set()) if e])
    
    # ===== 组装计量学信息 =====
    stoich_data: dict[str, dict] = {}
    all_reaction_nodes = set(reaction_node_to_id.keys())
    all_reaction_nodes.update(reaction_left_bnodes.keys())
    all_reaction_nodes.update(reaction_right_bnodes.keys())
    
    for rxn_node in all_reaction_nodes:
        mnxr = reaction_node_to_id.get(rxn_node)
        if mnxr is None:
            continue
        if target_mnxr_ids is not None and mnxr not in target_mnxr_ids:
            continue
        
        if mnxr not in stoich_data:
            stoich_data[mnxr] = {"reaction_nodes": {}}
        
        rxn_key = _last_fragment(str(rxn_node))
        if rxn_key not in stoich_data[mnxr]["reaction_nodes"]:
            stoich_data[mnxr]["reaction_nodes"][rxn_key] = {"left": [], "right": []}
        
        bucket = stoich_data[mnxr]["reaction_nodes"][rxn_key]
        
        # 处理 left 空白节点
        for bnode in reaction_left_bnodes.get(rxn_node, set()):
            chems = bnode_to_chems.get(bnode, set())
            comps = bnode_to_comps.get(bnode, set()) or {''}
            coef = bnode_to_coef.get(bnode, 1.0)
            for ch in chems:
                for cp in comps:
                    bucket["left"].append({"chem": ch, "comp": cp, "coef": float(coef)})
        
        # 处理 right 空白节点
        for bnode in reaction_right_bnodes.get(rxn_node, set()):
            chems = bnode_to_chems.get(bnode, set())
            comps = bnode_to_comps.get(bnode, set()) or {''}
            coef = bnode_to_coef.get(bnode, 1.0)
            for ch in chems:
                for cp in comps:
                    bucket["right"].append({"chem": ch, "comp": cp, "coef": float(coef)})
    
    # ===== 识别缺失信息的 MNXR =====
    missing_mnxr = []
    for mnxr_id in all_mnxr_ids:
        entry = annotations.get(mnxr_id)
        if entry is None or (
            entry["name"] is None 
            and not entry["xrefs"] 
            and not entry["ec"] 
            and not entry["raw_xrefs"]
        ):
            missing_mnxr.append(mnxr_id)
    
    return annotations, stoich_data, sorted(missing_mnxr)


def _extract_stoich_from_mnxmnet(
    file_path: str,
    target_mnxr_ids: set[str] | None = None,
) -> dict[str, dict]:
    """
    从 MNXmnet.ttl 提取反应的计量学信息（left/right + chem/comp/coef）
    """
    # 用于计量学信息的数据结构（类似 dump_mnxr_stoich.py）
    reaction_node_to_id: dict[URIRef, str] = {}
    reaction_left_nodes: dict[URIRef, set[URIRef]] = {}
    reaction_right_nodes: dict[URIRef, set[URIRef]] = {}
    reaction_side_nodes: dict[URIRef, set[URIRef]] = {}
    
    side_node_to_chems: dict[URIRef, set[str]] = {}
    side_node_to_comps: dict[URIRef, set[str]] = {}
    side_node_to_coef: dict[URIRef, float] = {}

    def on_triple(triple: tuple[object, object, object]) -> None:
        s, p, o = triple
        if not isinstance(s, URIRef):
            return
        
        su = str(s)
        pl = _pred_local_name(p)
        
        # ===== 处理反应的 mnxr 标记 =====
        if any(su.startswith(pref) for pref in REAC_SUBJ_PREFIXES):
            if pl == "mnxr":
                mnxr_id = _last_fragment(str(o))
                if mnxr_id.upper().startswith("MNXR"):
                    if target_mnxr_ids is None or mnxr_id in target_mnxr_ids:
                        reaction_node_to_id[s] = mnxr_id
        
        # ===== 处理反应的 left/right/side 节点 =====
        if pl == "left" and isinstance(o, URIRef):
            reaction_left_nodes.setdefault(s, set()).add(o)
        elif pl == "right" and isinstance(o, URIRef):
            reaction_right_nodes.setdefault(s, set()).add(o)
        elif pl == "side" and isinstance(o, URIRef):
            reaction_side_nodes.setdefault(s, set()).add(o)
        
        # ===== 处理 side 节点的属性（chem/comp/coef） =====
        if pl == "chem" and isinstance(o, URIRef):
            ou = str(o)
            if any(ou.startswith(ns) for ns in CHEM_NS_SET):
                side_node_to_chems.setdefault(s, set()).add(_last_fragment(ou))
        elif pl == "comp" and isinstance(o, URIRef):
            ou = str(o)
            if any(ou.startswith(ns) for ns in COMP_NS_SET):
                side_node_to_comps.setdefault(s, set()).add(_last_fragment(ou))
        elif pl == "coef":
            try:
                side_node_to_coef[s] = float(o)
            except Exception:
                pass

    # 解析（流式）
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
    
    # ===== 组装计量学信息 =====
    stoich_data: dict[str, dict] = {}
    all_reaction_nodes = set()
    all_reaction_nodes.update(reaction_left_nodes.keys())
    all_reaction_nodes.update(reaction_right_nodes.keys())
    all_reaction_nodes.update(reaction_side_nodes.keys())
    
    for rxn_node in all_reaction_nodes:
        mnxr = reaction_node_to_id.get(rxn_node)
        if mnxr is None:
            continue
        if target_mnxr_ids is not None and mnxr not in target_mnxr_ids:
            continue
        
        if mnxr not in stoich_data:
            stoich_data[mnxr] = {"reaction_nodes": {}}
        
        rxn_key = _last_fragment(str(rxn_node))
        if rxn_key not in stoich_data[mnxr]["reaction_nodes"]:
            stoich_data[mnxr]["reaction_nodes"][rxn_key] = {"left": [], "right": []}
        
        bucket = stoich_data[mnxr]["reaction_nodes"][rxn_key]
        
        # 处理 left
        for side_node in reaction_left_nodes.get(rxn_node, set()):
            chems = side_node_to_chems.get(side_node, set())
            comps = side_node_to_comps.get(side_node, set()) or {''}
            coef = side_node_to_coef.get(side_node, 1.0)
            for ch in chems:
                for cp in comps:
                    bucket["left"].append({"chem": ch, "comp": cp, "coef": float(coef)})
        
        # 处理 right
        for side_node in reaction_right_nodes.get(rxn_node, set()):
            chems = side_node_to_chems.get(side_node, set())
            comps = side_node_to_comps.get(side_node, set()) or {''}
            coef = side_node_to_coef.get(side_node, 1.0)
            for ch in chems:
                for cp in comps:
                    bucket["right"].append({"chem": ch, "comp": cp, "coef": float(coef)})
    
    return stoich_data


def sample_mnxref_content(file_path: str, max_reactions: int = 50) -> None:
    """
    采样 MNXref.ttl 的内容，输出少量反应的详细信息供检查
    """
    print(f"=== 正在采样 {file_path} 的前 {max_reactions} 个 MNXR 反应 ===\n")
    
    # 收集所有三元组（包括空白节点）
    mnxr_triples = {}  # mnxr_id -> [(pred, obj), ...]
    bnode_triples = {}  # bnode -> [(pred, obj), ...]
    bnodes_to_track = set()  # 需要追踪的空白节点
    count = 0
    stop_collecting_reactions = False
    
    def on_triple(triple: tuple[object, object, object]) -> None:
        nonlocal count, stop_collecting_reactions
        s, p, o = triple
        
        # 收集 MNXR 反应的三元组
        if isinstance(s, URIRef) and not stop_collecting_reactions:
            su = str(s)
            if any(su.startswith(pref) for pref in REAC_SUBJ_PREFIXES):
                mnxr_id = _last_fragment(su)
                if not mnxr_id.upper().startswith("MNXR"):
                    return
                
                if mnxr_id not in mnxr_triples:
                    mnxr_triples[mnxr_id] = []
                    count += 1
                    if count > max_reactions:
                        stop_collecting_reactions = True
                
                pred_name = _pred_local_name(p)
                obj_str = _last_fragment(str(o)) if isinstance(o, URIRef) else str(o)
                # 保留原始对象用于 BNode 匹配
                mnxr_triples[mnxr_id].append((pred_name, obj_str, o, type(o).__name__))
                
                # 记录需要追踪的空白节点
                if isinstance(o, BNode):
                    bnodes_to_track.add(o)
        
        # 收集空白节点的三元组（用于展开 left/right/side）
        if isinstance(s, BNode):
            pred_name = _pred_local_name(p)
            obj_str = _last_fragment(str(o)) if isinstance(o, URIRef) else str(o)
            # 使用 BNode 对象本身作为 key，而不是字符串
            if s not in bnode_triples:
                bnode_triples[s] = []
            bnode_triples[s].append((pred_name, obj_str))
    
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
    try:
        parser.parse(source, graph_like)
    except StopParsing:
        pass
    
    # 输出样本
    print(f"共采样到 {len(mnxr_triples)} 个 MNXR 反应，{len(bnode_triples)} 个空白节点")
    print(f"空白节点样例: {list(bnode_triples.keys())[:5]}\n")
    
    for idx, (mnxr_id, triples) in enumerate(list(mnxr_triples.items())[:3], 1):
        print(f"{'='*70}")
        print(f"反应 {idx}: {mnxr_id}")
        print(f"{'='*70}")
        
        # 按类别分组显示
        names = [obj for pred, obj, _, _ in triples if pred in {"name", "label"}]
        xrefs = [obj for pred, obj, _, _ in triples if "xref" in pred.lower() or pred in {"is", "seealso"}]
        ecs = [obj for pred, obj, _, _ in triples if "ec" in pred.lower()]
        # 保留原始对象用于匹配 BNode
        left_data = [(obj_str, orig_o, typ) for pred, obj_str, orig_o, typ in triples if pred == "left"]
        right_data = [(obj_str, orig_o, typ) for pred, obj_str, orig_o, typ in triples if pred == "right"]
        
        if names:
            print(f"  名称: {names[0][:60]}")
        if xrefs:
            print(f"  交叉引用: {', '.join(xrefs[:5])}")
        if ecs:
            print(f"  EC编号: {', '.join(ecs)}")
        
        print(f"\n  左侧节点 ({len(left_data)}):")
        for obj_str, orig_o, typ in left_data[:2]:  # 只显示前2个
            print(f"    - 节点类型: {typ}, 值: {obj_str[:40]}")
            if orig_o in bnode_triples:
                bn_info = bnode_triples[orig_o]
                chem = [o for p, o in bn_info if p == "chem"]
                comp = [o for p, o in bn_info if p == "comp"]
                coef = [o for p, o in bn_info if p == "coef"]
                if chem:
                    print(f"      → chem: {chem[0]}")
                if comp:
                    print(f"      → comp: {comp[0]}")
                if coef:
                    print(f"      → coef: {coef[0]}")
        
        print(f"\n  右侧节点 ({len(right_data)}):")
        for obj_str, orig_o, typ in right_data[:2]:  # 只显示前2个
            print(f"    - 节点类型: {typ}, 值: {obj_str[:40]}")
            if orig_o in bnode_triples:
                bn_info = bnode_triples[orig_o]
                chem = [o for p, o in bn_info if p == "chem"]
                comp = [o for p, o in bn_info if p == "comp"]
                coef = [o for p, o in bn_info if p == "coef"]
                if chem:
                    print(f"      → chem: {chem[0]}")
                if comp:
                    print(f"      → comp: {comp[0]}")
                if coef:
                    print(f"      → coef: {coef[0]}")
        print()
    
    # 统计分析
    total_with_stoich = sum(1 for triples in mnxr_triples.values()
                           if any(pred in {"left", "right", "side"} 
                                 for pred, _, _, _ in triples))
    print(f"\n{'='*70}")
    print(f"统计信息")
    print(f"{'='*70}")
    print(f"样本中有 left/right 的反应: {total_with_stoich}/{len(mnxr_triples)}")
    print(f"比例: {total_with_stoich/len(mnxr_triples)*100:.1f}%")
    print(f"\n✅ MNXref.ttl 包含计量学信息！")
    print(f"结构：反应 -> left/right 空白节点 -> chem/comp/coef")


if __name__ == "__main__":
    import csv
    import os
    
    mnxref_path = "/home/pickleopear/AlphaGEM/data_available/MNXref.ttl"
    mnxmnet_path = "/home/pickleopear/AlphaGEM/data_available/MNXmnet.ttl"
    out_dir = "/home/pickleopear/AlphaGEM/data_available"
    
    # 选择模式：
    MODE = "full"  # "sample" 或 "full"
    
    if MODE == "sample":
        # 采样模式：查看少量反应的详细信息
        print("正在采样 MNXref.ttl 内容...\n")
        sample_mnxref_content(mnxref_path, max_reactions=50)
        
        print("\n" + "="*60)
        print("请根据上述样本决定：")
        print("1. 如果 MNXref.ttl 包含计量学信息 -> 仅使用 MNXref.ttl")
        print("2. 如果 MNXref.ttl 不包含计量学信息 -> 同时使用 MNXref.ttl + MNXmnet.ttl")
        print("="*60)
    
    elif MODE == "full":
        # 完整模式：提取所有反应（注释 + 计量学信息）
        print("="*70)
        print("开始提取所有 MNXR 反应（注释 + 计量学信息）")
        print("="*70)
        
        # 使用组合函数同时提取注释和计量学信息
        ann, stoich, missing = read_mnxref_reacs_with_stoich(mnxref_path, mnxmnet_path)
        
        print(f"\n提取了 {len(ann)} 个反应的注释信息")
        print(f"提取了 {len(stoich)} 个反应的计量学信息")
        print(f"发现 {len(missing)} 个缺失详细信息的 MNXR")
        
        # 保存结果
        os.makedirs(out_dir, exist_ok=True)
        
        # 1. 完整注释（含 name/xrefs/ec/raw_xrefs）
        with open(f"{out_dir}/mnxref_reac_ann.pkl", "wb") as f:
        pickle.dump(ann, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"\n✓ 已保存注释到 {out_dir}/mnxref_reac_ann.pkl")
        
        # 2. 计量学信息
        with open(f"{out_dir}/mnxref_reac_stoich.pkl", "wb") as f:
            pickle.dump(stoich, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"✓ 已保存计量学信息到 {out_dir}/mnxref_reac_stoich.pkl")
        
        # 3. 仅 xrefs（向后兼容）
        with open(f"{out_dir}/mnxref_reac_xrefs.pkl", "wb") as f:
        pickle.dump({k: v.get("xrefs", {}) for k, v in ann.items()}, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        # 4. 原始 xref 三元组（可选）
        with open(f"{out_dir}/mnxref_reac_raw_xrefs.pkl", "wb") as f:
        pickle.dump({k: v.get("raw_xrefs", []) for k, v in ann.items()}, f, protocol=pickle.HIGHEST_PROTOCOL)

        # 5. 输出缺失注释信息的 MNXR 到 CSV
        if missing:
            missing_csv = f"{out_dir}/mnxref_reac_missing_annotation.csv"
            with open(missing_csv, "w", newline="") as fcsv:
                writer = csv.writer(fcsv)
                writer.writerow(["mnxr", "reason"])
                for mnxr in missing:
                    writer.writerow([mnxr, "no_annotation"])
            print(f"✓ 已输出缺失注释信息的 MNXR 到 {missing_csv}")
        
        # 6. 输出缺失计量学信息的 MNXR 到 CSV
        missing_stoich = []
        for mnxr in ann.keys():
            if mnxr not in stoich:
                missing_stoich.append(mnxr)
        
        if missing_stoich:
            missing_stoich_csv = f"{out_dir}/mnxref_reac_missing_stoich.csv"
            with open(missing_stoich_csv, "w", newline="") as fcsv:
                writer = csv.DictWriter(fcsv, fieldnames=["mnxr", "has_name", "has_xrefs", "has_ec"])
                writer.writeheader()
                for mnxr in missing_stoich:
                    entry = ann[mnxr]
                    writer.writerow({
                        "mnxr": mnxr,
                        "has_name": "yes" if entry.get("name") else "no",
                        "has_xrefs": "yes" if entry.get("xrefs") else "no",
                        "has_ec": "yes" if entry.get("ec") else "no"
                    })
            print(f"✓ 已输出缺失计量学信息的 MNXR 到 {missing_stoich_csv} ({len(missing_stoich)} 个)")
        
        print("\n" + "="*70)
        print("统计摘要")
        print("="*70)
        print(f"总反应数（有注释）: {len(ann)}")
        print(f"有计量学信息: {len(stoich)}")
        print(f"缺少注释信息: {len(missing)}")
        print(f"缺少计量学信息: {len(missing_stoich)}")
        print("="*70)
        print("完成！所有文件已保存")
        print("="*70)



