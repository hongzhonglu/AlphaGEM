from rdflib import Namespace, URIRef
import pickle
try:
    from rdflib.plugins.parsers.ttl import TurtleParser
except Exception:  # 兼容旧版本 rdflib
    from rdflib.plugins.parsers.notation3 import TurtleParser
from rdflib.parser import create_input_source


SCHEMA = Namespace("https://rdf.metanetx.org/schema/")
REAC = Namespace("https://rdf.metanetx.org/reac/")


class StopParsing(Exception):
    pass


def _extract_last_fragment(uri: URIRef) -> str:
    s = str(uri)
    if "#" in s:
        return s.rsplit("#", 1)[-1]
    return s.rstrip("/").rsplit("/", 1)[-1]


def extract_reaction_sides_from_turtle(
    file_path: str,
    target_mnxr_ids: set[str],
) -> dict[str, dict[str, list[str]]]:
    """
    流式解析 Turtle，不将整图载入内存；仅跟踪目标 MNXR 反应的左右底物。

    返回: { MNXRxx: {"left": [chem_uri...], "right": [chem_uri...] }, ... }
    """
    target = set(target_mnxr_ids)
    results: dict[str, dict[str, list[str]]] = {}

    # 仅追踪与目标反应相关的节点，降低内存占用
    reaction_node_to_id: dict[object, str] = {}
    reaction_left_nodes: dict[object, set[object]] = {}
    reaction_right_nodes: dict[object, set[object]] = {}
    side_node_to_chems: dict[object, list[str]] = {}

    def try_finalize(reaction_node: object) -> None:
        mnxr_id = reaction_node_to_id.get(reaction_node)
        if not mnxr_id or mnxr_id in results:
            return
        left_nodes = reaction_left_nodes.get(reaction_node)
        right_nodes = reaction_right_nodes.get(reaction_node)
        if not left_nodes or not right_nodes:
            return
        left_chems: list[str] = []
        right_chems: list[str] = []
        for ln in left_nodes:
            if ln in side_node_to_chems:
                left_chems.extend(side_node_to_chems[ln])
        for rn in right_nodes:
            if rn in side_node_to_chems:
                right_chems.extend(side_node_to_chems[rn])
        if left_chems and right_chems:
            results[mnxr_id] = {"left": left_chems, "right": right_chems}
            # 找齐所有目标后可提前终止解析
            if set(results.keys()) >= target:
                raise StopParsing()

    def on_triple(triple: tuple[object, object, object]) -> None:
        s, p, o = triple

        # 反应节点标识: ?rxn schema:mnxr <.../reac/MNXRxx>
        if p == SCHEMA.mnxr and isinstance(o, URIRef):
            mnxr_id = _extract_last_fragment(o)
            if mnxr_id in target:
                reaction_node_to_id[s] = mnxr_id
                try_finalize(s)
            return

        # 左右边：不依赖三元组出现顺序，先收集，若稍后识别为目标反应会被利用
        if p == SCHEMA.left:
            reaction_left_nodes.setdefault(s, set()).add(o)
            if s in reaction_node_to_id:
                try_finalize(s)
            return
        if p == SCHEMA.right:
            reaction_right_nodes.setdefault(s, set()).add(o)
            if s in reaction_node_to_id:
                try_finalize(s)
            return

        # 侧节点到化合物: ?side schema:chem <...chem/...>（顺序无关，全部暂存）
        if p == SCHEMA.chem and isinstance(o, URIRef):
            side_node_to_chems.setdefault(s, []).append(str(o))
            for rxn_node in list(reaction_node_to_id.keys()):
                try_finalize(rxn_node)

    parser = TurtleParser()
    try:
        source = create_input_source(location=file_path)
        # 最小 Graph-like：具备 store.context_aware、bind、namespace_manager、add、addN
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
            # 解析器可能写入四元组 (s,p,o,ctx)
            def addN(self, quads):
                for s, p, o, _ in quads:
                    on_triple((s, p, o))
            # 解析器需要：将相对 IRI 变为绝对 IRI
            def absolutize(self, uri, defrag=True):
                return urljoin(self._base, str(uri))
            # 某些分支会尝试替换 BNode；此处返回 None 表示不替换
            def skolem_genid(self, bnode):
                return None

        graph_like = _GraphLike()
        parser.parse(source, graph_like)
    except StopParsing:
        pass

    return results


def read_ttl_file(file_path: str, target_ids: list[str] | None = None):
    """
    示例：按需读取大体量 Turtle 文件，提取给定 MNXR 反应的左右底物。
    """
    if target_ids is None or len(target_ids) == 0:
        print("请提供至少一个 MNXR 反应 ID，例如: ['MNXR01']")
        return
    res = extract_reaction_sides_from_turtle(file_path, set(target_ids))
    for mnxr, sides in res.items():
        left_str = ", ".join(sides.get("left", []))
        right_str = ", ".join(sides.get("right", []))
        print(f"{mnxr}: LEFT[{left_str}] -> RIGHT[{right_str}]")


# ------------------------
# 全量解析所有 MNXR 反应
# ------------------------
def extract_all_mnxr_from_turtle(
    file_path: str,
) -> dict[str, dict[str, dict[str, list[str]]]]:
    """
    流式解析 Turtle，提取所有以 MNXR 开头的反应的左右底物与所属 comp。

    返回结构（按 MNXR 聚合，内部再按具体反应节点区分）：
    {
      MNXRxx: {
        REACTION_KEY: {
          "left": ["chem_uri@comp_uri", ...],
          "right": ["chem_uri@comp_uri", ...],
        },
        ...
      },
      ...
    }
    """
    # 关系缓存
    reaction_node_to_id: dict[object, str] = {}
    reaction_left_nodes: dict[object, set[object]] = {}
    reaction_right_nodes: dict[object, set[object]] = {}
    side_node_to_chems: dict[object, set[str]] = {}
    side_node_to_comps: dict[object, set[str]] = {}
    # 反向索引：side 节点 -> 关联的 reaction 节点
    side_node_to_reactions: dict[object, set[object]] = {}

    # 结果以 set 聚合，结尾再转 list
    # 第一层 MNXR -> 第二层 反应节点唯一键 -> {left/right: set("chem@comp")}
    accum: dict[str, dict[str, dict[str, set[str]]]] = {}

    def _node_key(n: object) -> str:
        return str(n)

    def ensure_bucket(mnxr_id: str, reaction_node: object) -> dict[str, set[str]]:
        if mnxr_id not in accum:
            accum[mnxr_id] = {}
        rkey = _node_key(reaction_node)
        if rkey not in accum[mnxr_id]:
            accum[mnxr_id][rkey] = {"left": set(), "right": set()}
        return accum[mnxr_id][rkey]

    def finalize_for_reaction(reaction_node: object) -> None:
        mnxr_id = reaction_node_to_id.get(reaction_node)
        if not mnxr_id:
            return
        left_nodes = reaction_left_nodes.get(reaction_node, set())
        right_nodes = reaction_right_nodes.get(reaction_node, set())
        if not left_nodes and not right_nodes:
            return
        bucket = ensure_bucket(mnxr_id, reaction_node)

        def add_pairs(side_node: object, side_key: str) -> None:
            chems = side_node_to_chems.get(side_node, set())
            comps = side_node_to_comps.get(side_node, set())
            # 仅取 URI 最后一段作为标识
            chem_ids = {_extract_last_fragment(c) for c in chems}
            comp_ids = {_extract_last_fragment(c) for c in comps} if comps else set()
            if comp_ids:
                for chem_id in chem_ids:
                    # 如果之前已添加过无 comp 的占位，先移除占位，避免重复
                    placeholder = f"{chem_id}@"
                    if placeholder in bucket[side_key]:
                        bucket[side_key].discard(placeholder)
                    for comp_id in comp_ids:
                        bucket[side_key].add(f"{chem_id}@{comp_id}")
            else:
                # 若未解析到 comp，仍保持 chem@ 格式（仅使用末段ID）
                for chem_id in chem_ids:
                    bucket[side_key].add(f"{chem_id}@")

        for ln in left_nodes:
            add_pairs(ln, "left")
        for rn in right_nodes:
            add_pairs(rn, "right")

    def on_triple(triple: tuple[object, object, object]) -> None:
        s, p, o = triple
        # 标记 MNXR 反应
        if p == SCHEMA.mnxr and isinstance(o, URIRef):
            mnxr_id = _extract_last_fragment(o)
            if mnxr_id.startswith("MNXR"):
                reaction_node_to_id[s] = mnxr_id
                finalize_for_reaction(s)
            return
        # 左/右边
        if p == SCHEMA.left:
            reaction_left_nodes.setdefault(s, set()).add(o)
            side_node_to_reactions.setdefault(o, set()).add(s)
            finalize_for_reaction(s)
            return
        if p == SCHEMA.right:
            reaction_right_nodes.setdefault(s, set()).add(o)
            side_node_to_reactions.setdefault(o, set()).add(s)
            finalize_for_reaction(s)
            return
        # 化合物
        if p == SCHEMA.chem and isinstance(o, URIRef):
            side_node_to_chems.setdefault(s, set()).add(str(o))
            # 回写所有关联反应
            for rxn in side_node_to_reactions.get(s, set()):
                finalize_for_reaction(rxn)
            return
        # compartment（若 schema 中存在）
        if getattr(SCHEMA, "comp", None) is not None and p == SCHEMA.comp and isinstance(o, URIRef):
            side_node_to_comps.setdefault(s, set()).add(str(o))
            for rxn in side_node_to_reactions.get(s, set()):
                finalize_for_reaction(rxn)
            return

    # 解析：与上面相同的最小 Graph-like 接收器
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
    parser.parse(source, graph_like)

    # 转 list 输出
    out: dict[str, dict[str, dict[str, list[str]]]] = {}
    for mnxr, per_reac in accum.items():
        out[mnxr] = {}
        for rkey, sides in per_reac.items():
            out[mnxr][rkey] = {
                "left": sorted(sides["left"]),
                "right": sorted(sides["right"]),
            }
    return out


def read_all_mnxr(file_path: str) -> dict[str, dict[str, dict[str, list[str]]]]:
    """
    便捷包装：返回所有 MNXR 反应，按反应节点区分，并以 chem@comp 形式存储左右。
    """
    return extract_all_mnxr_from_turtle(file_path)

if __name__ == "__main__":
    res = read_all_mnxr("/data_available/MNXmnet.ttl")
    # 保存为 pickle（仅存末段ID的 chem@comp）
    with open("/data_available/mnxr_left_right.pkl", "wb") as f:
        pickle.dump(res, f, protocol=pickle.HIGHEST_PROTOCOL)
    # 默认示例：仅解析与 MNXR01 相关的三元组，避免整文件载入
    # read_ttl_file("./data_available/MNXmnet.ttl", ["MNXR02", "MNXR03", "MNXR01", "MNXR152875"])