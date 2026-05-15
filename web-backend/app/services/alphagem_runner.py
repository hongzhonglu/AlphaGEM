from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from app.core.config import settings


def build_command(job: dict) -> list[str]:
    cmd = [
        "python",
        "AlphaGEM.py",
        "--name",
        job["internal_name"],
        "--mode",
        job["mode"],
        "--refname",
        job["refname"],
        "--fasta",
        job["fasta_path"],
        "--grothmedium",
        job["growth_medium"],
    ]

    if job.get("maplist_path"):
        cmd.extend(["--maplist", job["maplist_path"]])
    if job.get("structure_path"):
        cmd.extend(["--structure", job["structure_path"]])
    if job.get("tmscore") is not None:
        cmd.extend(["--TMscore", str(job["tmscore"])])
    if job.get("coverage") is not None:
        cmd.extend(["--coverage", str(job["coverage"])])
    if job.get("plddt") is not None:
        cmd.extend(["--pLDDT", str(job["plddt"])])

    return cmd


def resolve_output_dir(job: dict) -> Path:
    return settings.repo_root / "working" / job["internal_name"]


def collect_artifacts(output_dir: Path) -> list[Path]:
    if not output_dir.exists():
        return []

    candidates: list[Path] = []
    for path in output_dir.iterdir():
        if path.is_file():
            candidates.append(path)
    return sorted(candidates)


def copy_result_snapshot(output_dir: Path, snapshot_dir: Path) -> None:
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    if not output_dir.exists():
        return
    for item in output_dir.iterdir():
        destination = snapshot_dir / item.name
        if item.is_dir():
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(item, destination)
        else:
            shutil.copy2(item, destination)


def run_subprocess(job: dict, log_file: Path) -> int:
    cmd = build_command(job)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        process = subprocess.Popen(
            cmd,
            cwd=settings.repo_root,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return process.wait()
