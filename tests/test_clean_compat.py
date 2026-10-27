"""Check CLEAN's nested subprocess without running models or using a GPU."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import use_clean


def test_nested_esm_extraction_keeps_weights_only_and_interpreter(tmp_path, monkeypatch):
    app = tmp_path / "clean"
    (app / "src/CLEAN").mkdir(parents=True)
    (app / "src/CLEAN/__init__.py").write_text("")
    (app / "src/CLEAN/utils.py").write_text(
        "def retrive_esm1b_embedding(name):\n"
        "    raise RuntimeError('unwrapped ESM subprocess')\n"
    )
    (app / "esm/scripts").mkdir(parents=True)
    torch.save({"args": argparse.Namespace(layers=33)}, app / "allowed.pt")
    torch.save(Path("not_allowlisted"), app / "forbidden.pt")
    (app / "esm/scripts/extract.py").write_text(
        "import argparse, json, pickle, sys, torch\n"
        "from pathlib import Path\n"
        "assert torch.load('allowed.pt')['args'].layers == 33\n"
        "try:\n"
        "    torch.load('forbidden.pt')\n"
        "except pickle.UnpicklingError:\n"
        "    pass\n"
        "else:\n"
        "    raise AssertionError('weights-only loading was disabled')\n"
        "Path('extraction.json').write_text(json.dumps({\n"
        "    'python': sys.executable, 'args': sys.argv[1:]}))\n"
    )
    (app / "CLEAN_infer_fasta.py").write_text(
        "import sys\n"
        "from pathlib import Path\n"
        "from CLEAN.utils import retrive_esm1b_embedding\n"
        "name = sys.argv[2]\n"
        "retrive_esm1b_embedding('inputs/' + name)\n"
        "Path('results/inputs').mkdir(parents=True)\n"
        "Path('results/inputs/' + name + '_maxsep.csv').write_text('result')\n"
    )
    (tmp_path / "working/candida").mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(use_clean, "clean_path", str(app))

    use_clean.clean("candida")

    assert Path.cwd() == tmp_path
    assert (tmp_path / "working/candida/candida_homoleft_maxsep.csv").read_text() == "result"
    extraction = json.loads((app / "extraction.json").read_text())
    assert extraction["python"] == sys.executable
    assert extraction["args"] == ["esm1b_t33_650M_UR50S",
                                   "data/inputs/candida_homoleft.fasta",
                                   "data/esm_data", "--include", "mean"]


def test_clean_failure_does_not_move_stale_results(tmp_path, monkeypatch):
    app = tmp_path / "clean"
    (app / "src/CLEAN").mkdir(parents=True)
    (app / "src/CLEAN/__init__.py").write_text("")
    (app / "src/CLEAN/utils.py").write_text("")
    (app / "CLEAN_infer_fasta.py").write_text("raise RuntimeError('inference failed')\n")
    (app / "results/inputs").mkdir(parents=True)
    stale = app / "results/inputs/candida_homoleft_maxsep.csv"
    stale.write_text("old result")
    (tmp_path / "working/candida").mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(use_clean, "clean_path", str(app))

    with pytest.raises(subprocess.CalledProcessError):
        use_clean.clean("candida")

    assert Path.cwd() == tmp_path
    assert stale.read_text() == "old result"
    assert not (tmp_path / "working/candida/candida_homoleft_maxsep.csv").exists()
