"""PLM/sequence evidence thresholds, propagation, and empty-screen regressions."""

import argparse
import json
from pathlib import Path
import sys
import subprocess
from unittest.mock import Mock

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "plmsearchtools/plmsearch/plmsearch_util")]
from src import cov_filter, ss_predict_choose
from src.plmsearch_pipeline import (
    PlmsearchSettings, add_plmsearch_arguments, run_homolog_search,
)


def parse(*arguments):
    parser = argparse.ArgumentParser()
    add_plmsearch_arguments(parser)
    return PlmsearchSettings.from_args(parser.parse_args(arguments))


def test_defaults_restore_historical_legacy_settings():
    assert parse() == PlmsearchSettings(0.8, 0.5, 20, 0.25, 60, "direct")


@pytest.mark.parametrize("option,value", [
    ("--plm-threshold", "nan"), ("--plm-coverage", "1.01"),
    ("--bbh-coverage", "-0.1"), ("--plm-identity", "101"),
    ("--bbh-identity", "inf"), ("--bbh-direction", "both"),
])
def test_invalid_thresholds_fail_before_pipeline(option, value):
    with pytest.raises(SystemExit):
        parse(option, value)


def test_pair_filter_requires_both_coverages_and_identity(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "working/test"
    work.mkdir(parents=True)
    pd.DataFrame([
        ["good", "target", .7, .6, 30.1],
        ["short_reference", "target", .4, .9, 80],
        ["short_target", "target", .9, .4, 80],
        ["identity_boundary", "target", .9, .9, 30],
        ["coverage_boundary", "target", .5, .9, 80],
    ], columns=["Gene1", "Gene2", "COV1", "COV2", "Identity"]).to_csv(
        work / "test_ss_0.8_blast_results.csv", index=False)
    result = cov_filter.cov_filter("test", "yeast", .8, .5, 30)
    assert list(result.Gene1) == ["good"]


def test_no_plm_candidates_still_produces_usable_empty_tables(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "working/test"
    work.mkdir(parents=True)
    (tmp_path / "data_available").mkdir()
    (work / "test_ss_predictor_filter1.0").write_text("")
    (work / "test.fasta").write_text(">sp|target|target\nAAAA\n")
    (tmp_path / "data_available/yeast.fasta").write_text(">sp|ref|ref\nAAAA\n")
    pd.DataFrame({"Entry": ["ref"], "Gene Names (ordered locus)": ["Y1"]}).to_excel(
        tmp_path / "data_available/yeast.xlsx", index=False)
    ss_predict_choose.ss_predict_choose("test", "yeast", 1.0)
    cov_filter.result_blast("test", "yeast", 1.0)
    result = cov_filter.cov_filter("test", "yeast", 1.0, .5, 20)
    assert result.empty
    assert {"Gene1", "Gene2", "COV1", "COV2", "Identity"} <= set(result.columns)


def test_plm_threshold_and_named_reference_mapping(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "working/test"
    work.mkdir(parents=True)
    (tmp_path / "data_available").mkdir()
    pd.DataFrame({"Gene Names (ordered locus)": ["Y1"], "Entry": ["ref"]}).to_excel(
        tmp_path / "data_available/yeast.xlsx", index=False)
    (work / "test_ss_predictor_filter0.8").write_text(
        "sp|target|target\tsp|ref|ref\t0.9\n"
        "sp|boundary|boundary\tsp|ref|ref\t0.8\n")
    ss_predict_choose.ss_predict_choose("test", "yeast", .8)
    assert pd.read_excel(work / "test_ss_homolog0.8.xlsx").iloc[:, -2:].values.tolist() == [["Y1", "target"]]


def test_bbh_recomputes_instead_of_reusing_coverage_only_cache(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "working/test"
    work.mkdir(parents=True)
    (tmp_path / "data_available").mkdir()
    pd.DataFrame({"Entry": ["ref"], "Gene Names (ordered locus)": ["Y1"]}).to_excel(
        tmp_path / "data_available/yeast.xlsx", index=False)
    (work / "bbh_parsed_0.25.csv").write_text("stale cache\n")
    get = Mock()
    monkeypatch.setattr(cov_filter, "make_blast_db", Mock())
    monkeypatch.setattr(cov_filter, "get_bbh", get)
    monkeypatch.setattr(cov_filter, "bbh_analysis", Mock(return_value=pd.DataFrame(
        [["t1", "ref", 70, "<=>"], ["t2", "ref", 80, "->"], ["t3", "ref", 60, "<=>"]],
        columns=["gene", "subject", "PID", "BBH"])))
    cov_filter.bbh("test", "yeast", .8, .25, 60.0, "bidirect")
    get.assert_called_once()
    result = pd.read_excel(work / "test_ss0.8_bbh_pid60.0_cov0.25_bidirect.xlsx")
    assert result[["Gene1", "Gene2"]].values.tolist() == [["Y1", "t1"]]


def test_blast_failure_stops_before_stale_output_can_be_parsed(monkeypatch):
    run = Mock(side_effect=subprocess.CalledProcessError(1, "blastp"))
    monkeypatch.setattr(cov_filter.subprocess, "run", run)
    with pytest.raises(subprocess.CalledProcessError):
        cov_filter.run_blastp("target.fa", "reference.fa", "working", "data_available")


@pytest.mark.parametrize("module_name", ["model_build1", "model_build2"])
def test_transfer_does_not_create_fake_genes_from_locus_prefixes(module_name):
    from importlib import import_module
    from cobra.core.gene import GPR
    module = import_module("src." + module_name)
    rule = module.str_chanege_for_gpr("YER060W or (YER060W-A and OTHER)", "YER060W", ["target1", "target2"])
    assert GPR.from_string(rule).genes == {"target1", "target2", "YER060W-A", "OTHER"}
    rule = module.str_chanege_for_gpr(rule, "YER060W-A", ["target3"])
    assert GPR.from_string(rule).genes == {"target1", "target2", "target3", "OTHER"}


@pytest.mark.parametrize("fail", [False, True])
def test_all_cli_settings_reach_stages_and_manifest(tmp_path, monkeypatch, fail):
    from src import generate_list, mkdir, plmsearch
    monkeypatch.chdir(tmp_path)
    settings = parse("--plm-threshold", ".85", "--plm-coverage", ".6", "--plm-identity", "35",
                     "--bbh-coverage", ".7", "--bbh-identity", "55", "--bbh-direction", "bidirect")
    work = tmp_path / "working/test"
    work.mkdir(parents=True)
    for path in [work / "test.fasta", tmp_path / "data_available/yeast.fasta",
                 tmp_path / "data_available/yeast.xlsx", work / "test_embedding.pkl",
                 tmp_path / "working/yeast_embedding.pkl", tmp_path / "plmsearchtools/model/plmsearch.sav"]:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"fixture")
    pd.DataFrame({0: ["Y1"], 1: ["target"]}).to_excel(work / "matrix_homologtest.xlsx")
    stages = []
    for module, name in [(generate_list, "generare_list_without_structure"), (plmsearch, "ss_predictor"),
                         (ss_predict_choose, "ss_predict_choose"), (cov_filter, "result_blast"),
                         (cov_filter, "cov_filter"), (cov_filter, "bbh"), (cov_filter, "merge"), (mkdir, "mvfile")]:
        mock = Mock(return_value=pd.DataFrame({"Gene1": ["Y1"]}))
        monkeypatch.setattr(module, name, mock)
        stages.append(mock)
    if fail:
        stages[4].side_effect = RuntimeError("screen failed")
        with pytest.raises(RuntimeError, match="screen failed"):
            run_homolog_search("test", "yeast", settings)
    else:
        run_homolog_search("test", "yeast", settings)
        stages[1].assert_called_once_with("yeast", "test", .85)
        stages[4].assert_called_once_with("test", "yeast", .85, .6, 35)
        stages[5].assert_called_once_with("test", "yeast", .85, .7, 55, "bidirect")
        for stage in stages[6:]:
            stage.assert_called_once_with("test", "yeast", .85, .6, 35, .7, 55, "bidirect")
    manifest = json.loads((work / "plmsearch_parameters.json").read_text())
    assert manifest["parameters"] == settings.__dict__
    assert manifest["status"] == ("failed" if fail else "completed")
    assert len(manifest["input_sha256"]) == 6
    if not fail:
        assert manifest["counts"]["merged_pairs"] == 1
