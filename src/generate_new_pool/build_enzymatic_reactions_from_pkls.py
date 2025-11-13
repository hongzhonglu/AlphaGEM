import os
import csv
import pickle
from typing import Dict, List, Set

import cobra
from cobra import Model, Reaction
from cobra.io import read_sbml_model, write_sbml_model


CHEM_SBML = "/home/pickleopear/AlphaGEM/data_available/chem.xml"
EC2MNXR_FULL = "/home/pickleopear/AlphaGEM/data_available/ec_to_mnxr_full.pkl"
EC2MNXR_GENERIC = "/home/pickleopear/AlphaGEM/data_available/ec_to_mnxr_generic.pkl"
REAC_XREFS_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxref_reac_xrefs.pkl"
REAC_ANN_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxref_reac_ann.pkl"
MNXR_STOICH_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxr_left_right.pkl"
MNXR_REVERSIBILITY_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxr_reversibility.pkl"
OUT_SBML = "/home/pickleopear/AlphaGEM/data_available/enzymatic_reactions.xml"
MISSING_MET_CSV = "/home/pickleopear/AlphaGEM/data_available/missing_metabolites.csv"
WATER_MET_ID = "WATER"
NO_EC_LIST = "/home/pickleopear/AlphaGEM/data_available/reactions_without_ec.txt"


DB_KEY_MAP: Dict[str, str] = {
    "KEGG": "kegg.reaction",
    "RHEA": "rhea",
    "METACYC": "metacyc.reaction",
    "REACTOME": "reactome",
    "BIGG": "bigg.reaction",
    "SEED": "seed.reaction",
}


def load_pickle(path: str) -> dict:
    if not os.path.exists(path):
        return {}
    with open(path, "rb") as f:
        return pickle.load(f)


def invert_ec_map(ec2mnxr: Dict[str, List[str]]) -> Dict[str, Set[str]]:
    reacs: Dict[str, Set[str]] = {}
    for ec, rxns in (ec2mnxr or {}).items():
        for r in rxns or []:
            reacs.setdefault(r, set()).add(ec)
    return reacs


def build_model_with_metabolites() -> Model:
    base = read_sbml_model(CHEM_SBML)
    model = Model("Enzymatic_Reactions")
    # 复制代谢物到新模型
    if base.metabolites:
        model.add_metabolites(list(base.metabolites))
    # 设定 compartments（统一为 c）
    model.compartments = {"c": "cytosol"}
    return model


def _get_or_create_water(model: Model):
    try:
        return model.metabolites.get_by_id(WATER_MET_ID)
    except KeyError:
        from cobra import Metabolite
        m = Metabolite(id=WATER_MET_ID, name="water", compartment="c")
        # 常用水的化学式与电荷
        m.formula = "H2O"
        m.charge = 0
        # 常见数据库标识
        m.annotation = {
            "chebi": ["CHEBI:15377"],
            "kegg.compound": ["C00001"],
            "pubchem.compound": ["962"],
            "hmdb": ["HMDB0002111"],
            "inchi": ["InChI=1S/H2O/h1H2"],
            "inchikey": ["XLYOFNOQVPJJNP-UHFFFAOYSA-N"],
        }
        model.add_metabolites([m])
        return m


def build_reaction(reac_id: str, xrefs: Dict[str, List[str]] | None, ecs: List[str] | None, name: str | None) -> Reaction:
    r = Reaction(id=reac_id, name=(name or reac_id))
    r.lower_bound = -1000.0
    r.upper_bound = 1000.0
    # 仅注释（无计量项，因当前仅基于 pkl，不解析左右底物）
    ann: Dict[str, List[str]] = {}
    for db_upper, ids in (xrefs or {}).items():
        key = DB_KEY_MAP.get(db_upper, db_upper.lower())
        if ids:
            ann[key] = list(ids)
    if ecs:
        ann["ec-code"] = sorted(list(set(ecs)))
    if ann:
        r.annotation = ann
    return r


def _parse_chem_comp(token: str) -> tuple[str, str | None]:
    # 形如 "MNXM123@MNXC3" 或仅 "MNXM123"
    if "@" in token:
        chem, comp = token.split("@", 1)
        chem = chem.strip()
        comp = comp.strip() or None
        return chem, comp
    return token.strip(), None


def _parse_item_to_coeff(item, default_coeff: float) -> tuple[str | None, float, str | None]:
    # item 可以是：
    # - 字符串："MNXM123@MNXC3" 或 "2.0|MNXM123@MNXC3" 或 "MNXM123@MNXC3|2"
    # - 字典：{"chem": "MNXM123@MNXC3", "coeff": 2.0} 或 {"met": ..., "stoich": ...}
    if isinstance(item, dict):
        chem_token = item.get("chem") or item.get("met") or item.get("metabolite") or item.get("id")
        # 支持多种命名，优先使用 coef（与 TTL 抽取保持一致）
        coeff = (
            item.get("coef")
            if item.get("coef") is not None
            else item.get("coeff") or item.get("stoich") or item.get("stoich_coeff") or item.get("n")
        )
        if chem_token is None:
            return None, default_coeff, None
        try:
            coeff_val = float(coeff) if coeff is not None else default_coeff
        except Exception:
            coeff_val = default_coeff
        chem_id, comp = _parse_chem_comp(str(chem_token))
        return chem_id, coeff_val, comp
    if isinstance(item, (list, tuple)) and len(item) >= 2:
        chem_token = str(item[0])
        chem_id, comp = _parse_chem_comp(chem_token)
        try:
            coeff_val = float(item[1])
        except Exception:
            coeff_val = default_coeff
        return chem_id, coeff_val, comp
    # 字符串：尝试用分隔符解析
    s = str(item)
    if "|" in s:
        parts = s.split("|")
        # 允许 "coef|chem" 或 "chem|coef"
        try:
            cval = float(parts[0])
            chem_id, comp = _parse_chem_comp(parts[1])
            return chem_id, cval, comp
        except Exception:
            pass
        try:
            cval = float(parts[-1])
            chem_id, comp = _parse_chem_comp("|".join(parts[:-1]))
            return chem_id, cval, comp
        except Exception:
            pass
    # 兜底：仅 chem，系数使用默认值
    chem_id, comp = _parse_chem_comp(s)
    return chem_id, default_coeff, comp


def _add_stoichiometry_from_pkl(
    r: Reaction,
    model: Model,
    side_items,
    is_left: bool,
    missing_rows: List[Dict[str, str]],
    mnxr: str,
    subkey: str | None,
    side_label: str,
) -> None:
    if not side_items:
        return
    sign = -1.0 if is_left else 1.0
    for item in side_items:
        chem_id, coeff, comp = _parse_item_to_coeff(item, default_coeff=1.0)
        if not chem_id:
            continue
        try:
            met = model.metabolites.get_by_id(chem_id)
        except KeyError:
            # 特殊处理 WATER：自动创建并加入模型
            if chem_id == WATER_MET_ID:
                met = _get_or_create_water(model)
            else:
                # 非 WATER：记录缺失并跳过
                missing_rows.append({
                    "reaction_id": r.id,
                    "mnxr": mnxr,
                    "subkey": str(subkey) if subkey is not None else "",
                    "side": side_label,
                    "chem_id": chem_id,
                    "comp": comp or "",
                    "coef": str(coeff),
                })
                continue
        r.add_metabolites({met: sign * float(coeff)})


def main():
    # 载入映射与注释
    ec_full = load_pickle(EC2MNXR_FULL)
    ec_generic = load_pickle(EC2MNXR_GENERIC)
    reac_xrefs = load_pickle(REAC_XREFS_PKL)  # MNXR -> {db:[ids]}
    mnxr_stoich = load_pickle(MNXR_STOICH_PKL)  # 结构：MNXR -> per_reaction_key -> {left:[...], right:[...]}
    mnxr_rev = load_pickle(MNXR_REVERSIBILITY_PKL)  # MNXR -> {direction: both/L2R/R2L, reversible: bool|None}
    reac_ann = load_pickle(REAC_ANN_PKL)  # MNXR -> {name, xrefs, ec}

    # 反转为 MNXR -> EC 列表
    rxn_to_ec: Dict[str, Set[str]] = invert_ec_map(ec_full)
    # 可选择加入 generic 的 EC（提升覆盖率）
    # gen_inv = invert_ec_map(ec_generic)
    # for r, ecs in gen_inv.items():
    #     rxn_to_ec.setdefault(r, set()).update(ecs)

    # 构建模型：先装入代谢物
    model = build_model_with_metabolites()
    missing_rows: List[Dict[str, str]] = []
    no_ec_rids: List[str] = []

    # 目标集合：带 EC 的 MNXR ∪ 具有计量信息的 MNXR（即使无 EC 也纳入）
    target_mnxr = set(rxn_to_ec.keys())
    if isinstance(mnxr_stoich, dict):
        target_mnxr |= set(mnxr_stoich.keys())

    reactions: List[Reaction] = []
    for mnxr in sorted(target_mnxr):
        ecs = list(rxn_to_ec.get(mnxr, set()))
        xref = reac_xrefs.get(mnxr) if isinstance(reac_xrefs, dict) else None
        rname = None
        if isinstance(reac_ann, dict):
            ra = reac_ann.get(mnxr)
            if isinstance(ra, dict):
                rname = ra.get("name")
        # 有可能一个 MNXR 对应多个具体反应节点
        per_reac = mnxr_stoich.get(mnxr, {}) if isinstance(mnxr_stoich, dict) else {}
        if isinstance(per_reac, dict) and per_reac:
            # 仅取第一个子反应（按键的字符串排序确定性选择）
            subkey, sides = next(iter(sorted(per_reac.items(), key=lambda kv: str(kv[0]))))
            # 使用原始 MNXR 作为反应 ID
            rid = mnxr
            r = build_reaction(rid, xref, list(ecs), rname)
            # 根据左右两侧补充计量
            left_items = (sides or {}).get("left")
            right_items = (sides or {}).get("right")
            _add_stoichiometry_from_pkl(
                r, model, left_items, is_left=True,
                missing_rows=missing_rows, mnxr=mnxr, subkey=subkey, side_label="left",
            )
            _add_stoichiometry_from_pkl(
                r, model, right_items, is_left=False,
                missing_rows=missing_rows, mnxr=mnxr, subkey=subkey, side_label="right",
            )
            # 根据可逆性设置边界
            rev = mnxr_rev.get(mnxr) if isinstance(mnxr_rev, dict) else None
            if isinstance(rev, dict):
                d = rev.get("direction")
                if d == "L2R":
                    r.lower_bound = 0.0
                elif d == "R2L":
                    r.upper_bound = 0.0
                elif d == "both":
                    r.lower_bound = -1000.0
                    r.upper_bound = 1000.0
            if not ecs:
                no_ec_rids.append(r.id)
            reactions.append(r)
        else:
            # 无详细分支，创建汇总反应但无计量（至少带注释）
            r = build_reaction(mnxr, xref, list(ecs), rname)
            rev = mnxr_rev.get(mnxr) if isinstance(mnxr_rev, dict) else None
            if isinstance(rev, dict):
                d = rev.get("direction")
                if d == "L2R":
                    r.lower_bound = 0.0
                elif d == "R2L":
                    r.upper_bound = 0.0
            if not ecs:
                no_ec_rids.append(r.id)
            reactions.append(r)

    if reactions:
        model.add_reactions(reactions)

    # 导出 SBML
    write_sbml_model(model, OUT_SBML)

    # 如有缺失代谢物，输出 CSV 提示
    if missing_rows:
        fieldnames = ["reaction_id", "mnxr", "subkey", "side", "chem_id", "comp", "coef"]
        with open(MISSING_MET_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(missing_rows)

    # 输出无 EC 的反应列表
    if no_ec_rids:
        with open(NO_EC_LIST, "w") as f:
            for rid in no_ec_rids:
                f.write(f"{rid}\n")
        print(f"无 EC 的反应数量: {len(no_ec_rids)} (已写入 {NO_EC_LIST})")
    else:
        print("无 EC 的反应数量: 0")


if __name__ == "__main__":
    main()


