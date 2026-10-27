"""Missing table rows are errors; documented proteins without loci are audited."""

import json
from pathlib import Path
import sys
from unittest.mock import Mock

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src'), str(ROOT / 'plmsearchtools/plmsearch/plmsearch_util')]
from src import cov_filter, model_build1, ss_predict_choose
from src.reference_mapping import load_reference_loci, map_reference_pairs


@pytest.fixture
def work(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'data_available').mkdir()
    work = tmp_path / 'working/test'
    work.mkdir(parents=True)
    pd.DataFrame({'Entry': ['valid', 'unassigned', 'spaces'],
                  'Gene Names (ordered locus)': ['Y1', None, '   ']}).to_excel(
        tmp_path / 'data_available/yeast.xlsx', index=False)
    return work


def test_blank_loci_are_excluded_with_pair_audit(work):
    frame = pd.DataFrame({'ref': ['valid', 'unassigned', 'spaces'], 'target': ['t1', 't2', 't3']})
    result = map_reference_pairs(frame, 'ref', 'target', 'yeast', work / 'audit.json')
    assert result.to_dict('records') == [{'ref': 'Y1', 'target': 't1'}]
    audit = json.loads((work / 'audit.json').read_text())
    assert audit['excluded_pair_count'] == 2
    assert audit['excluded_pairs'] == [
        {'reference_accession': 'unassigned', 'target_gene': 't2'},
        {'reference_accession': 'spaces', 'target_gene': 't3'},
    ]


def test_missing_accession_still_fails(work):
    frame = pd.DataFrame({'ref': ['absent'], 'target': ['t1']})
    with pytest.raises(ValueError, match='missing ordered locus mapping'):
        map_reference_pairs(frame, 'ref', 'target', 'yeast', work / 'audit.json')


@pytest.mark.parametrize('problem', ['missing_table', 'missing_column', 'duplicate_entry', 'blank_entry'])
def test_incomplete_reference_resources_fail_explicitly(work, problem):
    path = Path('data_available/yeast.xlsx')
    if problem == 'missing_table':
        path.unlink()
    elif problem == 'missing_column':
        pd.DataFrame({'Entry': ['valid']}).to_excel(path, index=False)
    else:
        pd.DataFrame({'Entry': ['valid', 'valid' if problem == 'duplicate_entry' else None],
                      'Gene Names (ordered locus)': ['Y1', 'Y2']}).to_excel(path, index=False)
    with pytest.raises((ValueError, FileNotFoundError)):
        load_reference_loci('yeast')


def test_plm_output_does_not_contain_nan_genes(work):
    (work / 'test_ss_predictor_filter0.8').write_text(
        'sp|t1|t1\tsp|valid|valid\t0.99\n'
        'sp|t2|t2\tsp|unassigned|unassigned\t0.99\n')
    ss_predict_choose.ss_predict_choose('test', 'yeast', .8)
    frame = pd.read_excel(work / 'test_ss_homolog0.8.xlsx').iloc[:, -2:]
    assert frame.values.tolist() == [['Y1', 't1']]
    assert json.loads((work / 'test_ss_homolog0.8_unmapped_reference.json').read_text())['excluded_pair_count'] == 1


def test_bbh_uses_same_no_locus_policy(work, monkeypatch):
    monkeypatch.setattr(cov_filter, 'make_blast_db', Mock())
    monkeypatch.setattr(cov_filter, 'get_bbh', Mock())
    monkeypatch.setattr(cov_filter, 'bbh_analysis', Mock(return_value=pd.DataFrame(
        [['t1', 'valid', 99, '<=>'], ['t2', 'unassigned', 99, '<=>']],
        columns=['gene', 'subject', 'PID', 'BBH'])))
    cov_filter.bbh('test', 'yeast', .8)
    frame = pd.read_excel(work / 'test_ss0.8_bbh_pid60.0_cov0.25_direct.xlsx')
    assert frame[['Gene1', 'Gene2']].values.tolist() == [['Y1', 't1']]
    assert json.loads((work / 'test_bbh_unmapped_reference.json').read_text())['excluded_pair_count'] == 1


def test_draft_export_removes_orphan_reference_genes(work, monkeypatch):
    import cobra

    model = cobra.Model('reference')
    model.solver = 'glpk'
    reaction = cobra.Reaction('R')
    reaction.gene_reaction_rule = 'Y1'
    model.add_reactions([reaction])
    pd.DataFrame({'reference': ['Y1'], 'target': ['target']}).to_excel(
        work / 'matrix_homologtest.xlsx')
    monkeypatch.setattr(cobra.io, 'read_sbml_model', lambda path: model)
    model_build1.modelbuild('yeast-GEM.xml', 'test')
    output = cobra.io.load_yaml_model(work / 'model1test.yml')
    assert {gene.id for gene in output.genes} == {'target'}
    assert output.reactions.R.gpr.genes == {'target'}
