"""Parameters and provenance for the PLM-supported homolog transfer route."""

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path


def fraction(value):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise argparse.ArgumentTypeError("must be a finite fraction between 0 and 1")
    return value


def percent(value):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 100:
        raise argparse.ArgumentTypeError("must be a finite percentage between 0 and 100")
    return value


def add_plmsearch_arguments(parser):
    group = parser.add_argument_group("PLMSearch homolog screening")
    group.add_argument("--plm-threshold", "--plm-score", type=fraction, default=0.8,
                       help="minimum PLM homolog score (strict >); default: 0.8")
    group.add_argument("--plm-coverage", type=fraction, default=0.5,
                       help="minimum BLAST coverage of BOTH PLM candidate proteins (strict >); default: 0.5")
    group.add_argument("--plm-identity", type=percent, default=20.0,
                       help="minimum BLAST identity for PLM candidates, in percent (strict >); default: 20")
    group.add_argument("--bbh-coverage", type=fraction, default=0.25,
                       help="minimum query coverage in each BLAST direction (>=); default: 0.25")
    group.add_argument("--bbh-identity", type=percent, default=60.0,
                       help="minimum BLAST best-hit identity, in percent (strict >); default: 60")
    group.add_argument("--bbh-direction", choices=["direct", "bidirect"], default="direct",
                       help="include one-way best hits or require reciprocal best hits; default: direct")


@dataclass(frozen=True)
class PlmsearchSettings:
    threshold: float = 0.8
    coverage: float = 0.5
    identity: float = 20.0
    bbh_coverage: float = 0.25
    bbh_identity: float = 60.0
    bbh_direction: str = "direct"

    @classmethod
    def from_args(cls, args):
        return cls(args.plm_threshold, args.plm_coverage, args.plm_identity,
                   args.bbh_coverage, args.bbh_identity, args.bbh_direction)


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_homolog_search(name, refname, settings):
    from . import cov_filter, generate_list, mkdir, plmsearch, ss_predict_choose
    import pandas as pd

    work = Path("working") / name
    manifest_path = work / "plmsearch_parameters.json"
    manifest = {
        "status": "running", "name": name, "reference": refname,
        "parameters": asdict(settings),
        "input_sha256": {
            str(path): file_sha256(path) for path in (
                work / f"{name}.fasta", Path(f"data_available/{refname}.fasta"),
                Path(f"data_available/{refname}.xlsx"), work / f"{name}_embedding.pkl",
                Path(f"working/{refname}_embedding.pkl"),
                Path("plmsearchtools/model/plmsearch.sav"),
            )
        },
        "screening": "union of PLM+pairwise BLAST and independent BLAST best-hit evidence",
    }

    def write_manifest():
        temporary = manifest_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(manifest, indent=2) + "\n")
        temporary.replace(manifest_path)

    write_manifest()
    score, cov, identity = settings.threshold, settings.coverage, settings.identity
    bcov, bid, direction = settings.bbh_coverage, settings.bbh_identity, settings.bbh_direction
    print("PLMSearch homolog parameters:", json.dumps(asdict(settings)), flush=True)
    try:
        generate_list.generare_list_without_structure(name)
        plmsearch.ss_predictor(refname, name, score)
        ss_predict_choose.ss_predict_choose(name, refname, score)
        cov_filter.result_blast(name, refname, score)
        plm_hits = cov_filter.cov_filter(name, refname, score, cov, identity)
        cov_filter.bbh(name, refname, score, bcov, bid, direction)
        cov_filter.merge(name, refname, score, cov, identity, bcov, bid, direction)
        mkdir.mvfile(name, refname, score, cov, identity, bcov, bid, direction)
        matrix_path = work / f"matrix_homolog{name}.xlsx"
        matrix = pd.read_excel(matrix_path).iloc[:, -2:].drop_duplicates()
        manifest.update(status="completed", counts={
            "plm_blast_pairs": len(plm_hits), "merged_pairs": len(matrix),
            "reference_genes": int(matrix.iloc[:, 0].nunique()),
            "target_genes": int(matrix.iloc[:, 1].nunique()),
        }, matrix_sha256=file_sha256(matrix_path))
    except Exception as error:
        manifest.update(status="failed", error=str(error))
        raise
    finally:
        write_manifest()
    return manifest
