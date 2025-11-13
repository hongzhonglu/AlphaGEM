import os
import sys
import pickle
from typing import Dict, List
sys.path.append('./src')
import cobra
from cobra import Model, Metabolite
from cobra.io import write_sbml_model


DEFAULT_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxref_chem_annotations.pkl"
DEFAULT_OUT = "/home/pickleopear/AlphaGEM/data_available/chem.xml"


DB_KEY_MAP: Dict[str, str] = {
    "CHEBI": "chebi",
    "KEGG": "kegg.compound",
    "HMDB": "hmdb",
    "PUBCHEM": "pubchem.compound",
    "METACYC": "metacyc",
    "BIGG": "bigg.metabolite",
    "SEED": "seed.compound",
    "INCHI": "inchi",
    "INCHIKEY": "inchikey",
    # 其余未映射的键将按小写原样写入
}


def build_metabolite(mnxm_id: str, info: Dict) -> Metabolite:
    # 使用描述性名称：优先 description，其次 name/label，不再回退为 id
    name = info.get("description") or info.get("name") or info.get("label")
    formula = info.get("formula") or None
    charge_val = info.get("charge")
    try:
        charge = int(charge_val) if charge_val is not None else None
    except Exception:
        charge = None

    m = Metabolite(id=mnxm_id, name=(name or ""), compartment="c")
    if formula:
        m.formula = formula
    if charge is not None:
        m.charge = charge

    # 注释：从 xrefs 整理到 MIRIAM 风格键
    xrefs: Dict[str, List[str]] = {}
    raw = info.get("xrefs") or {}
    for db_upper, ids in raw.items():
        key = DB_KEY_MAP.get(db_upper, db_upper.lower())
        if not ids:
            continue
        xrefs[key] = list(ids)
    # 始终保留 MNX 映射，避免 name 不含 id 时难以回溯
    xrefs.setdefault("metanetx.chemical", [mnxm_id])
    if xrefs:
        m.annotation = xrefs

    return m


def export_all_chems_to_sbml(pkl_path: str = DEFAULT_PKL, out_path: str = DEFAULT_OUT) -> None:
    with open(pkl_path, "rb") as f:
        ann = pickle.load(f)

    model = Model("MNXref_CHEM_only")
    model.compartments = {"c": "cytosol"}

    mets: List[Metabolite] = []
    for mnxm_id, info in ann.items():
        if not mnxm_id.startswith("MNXM"):
            continue
        # 排除 smiles，仅填充 name/formula/charge 与 xrefs
        m = build_metabolite(mnxm_id, info)
        mets.append(m)

    if mets:
        model.add_metabolites(mets)

    write_sbml_model(model, out_path)


if __name__ == "__main__":
    export_all_chems_to_sbml()


