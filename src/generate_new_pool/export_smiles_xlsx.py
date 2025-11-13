import os
import pickle
from typing import Dict, Any, List

import pandas as pd


DEFAULT_PKL = "/home/pickleopear/AlphaGEM/data_available/mnxref_chem_annotations.pkl"
DEFAULT_CSV = "/home/pickleopear/AlphaGEM/data_available/chem_smiles.csv"


def export_smiles_csv(pkl_path: str = DEFAULT_PKL, out_path: str = DEFAULT_CSV) -> None:
    with open(pkl_path, "rb") as f:
        ann: Dict[str, Dict[str, Any]] = pickle.load(f)

    rows: List[Dict[str, Any]] = []
    for mnxm, info in ann.items():
        if not isinstance(mnxm, str) or not mnxm.startswith("MNXM"):
            continue
        smiles = (info or {}).get("smiles")
        rows.append({
            "mnxm": mnxm,
            "smiles": smiles if smiles is not None else "",
        })

    df = pd.DataFrame(rows).sort_values(by=["mnxm"]).reset_index(drop=True)
    # 确保目录存在
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)


if __name__ == "__main__":
    export_smiles_csv()


