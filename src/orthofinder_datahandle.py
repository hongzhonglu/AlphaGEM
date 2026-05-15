import os
import re
import shutil
import subprocess
from itertools import product
from pathlib import Path

import pandas as pd


MAX_GROUP_MEMBERS = 200
ORTHOFINDER_ID_PATTERN = re.compile(r"\|([^|]+)\|")


def _extract_orthofinder_ids(cell):
    if not isinstance(cell, str) or not cell:
        return []
    return ORTHOFINDER_ID_PATTERN.findall(cell)


def _latest_orthogroups_tsv(orthofinder_dir):
    result_dirs = [item for item in orthofinder_dir.iterdir() if item.is_dir()]
    if not result_dirs:
        raise FileNotFoundError(f"No OrthoFinder result directory found in {orthofinder_dir}")
    latest_dir = max(result_dirs, key=lambda item: item.stat().st_mtime)
    return latest_dir / "Orthogroups" / "Orthogroups.tsv"


def datahandel(name="", refname=""):
    path = Path.cwd()
    species_dir = path / "orth" / "data" / name

    if species_dir.exists():
        shutil.rmtree(species_dir)
    species_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy(path / "data_available" / f"{refname}.fasta", species_dir / f"z-{refname}.fasta")
    shutil.copy(path / "working" / name / f"{name}.fasta", species_dir / f"a-{name}.fasta")

    subprocess.run(
        ["python",str(path / "tools" / "OrthoFinder" / "orthofinder.py"), "-f", str(species_dir), "-og"],
        check=True,
    )

    pathresult = _latest_orthogroups_tsv(species_dir / "OrthoFinder")
    orth = pd.read_csv(pathresult, sep="\t", dtype=str).fillna("")

    pair_rows = []
    for _, row in orth.iterrows():
        target_ids = _extract_orthofinder_ids(row.iloc[1])[:MAX_GROUP_MEMBERS]
        ref_ids = _extract_orthofinder_ids(row.iloc[2])[:MAX_GROUP_MEMBERS]
        if not target_ids or not ref_ids:
            continue
        pair_rows.extend(product(ref_ids, target_ids))

    juzhen1 = pd.DataFrame(pair_rows)

    ref = pd.read_excel(path / "data_available" / f"{refname}.xlsx", dtype=str)
    ref_map = dict(zip(ref.iloc[:, 0], ref.iloc[:, 2]))
    if not juzhen1.empty:
        juzhen1.iloc[:, 0] = juzhen1.iloc[:, 0].map(lambda value: ref_map.get(value, value))

    juzhen1.to_excel(path / "working" / name / f"matrix_orthofinder{name}.xlsx")
