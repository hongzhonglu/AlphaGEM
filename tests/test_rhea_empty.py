"""The actual Rhea XLSX producer can emit a workbook with no data columns."""

from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src'), str(ROOT / 'plmsearchtools/plmsearch/plmsearch_util')]
from src.use_rhea import process_rhea


@pytest.mark.parametrize('kind', ['no_columns', 'header_only', 'hits'])
def test_rhea_producer_workbooks_are_read_without_dropping_real_hits(tmp_path, monkeypatch, kind):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / 'working/test'
    work.mkdir(parents=True)
    (tmp_path / 'rhea').mkdir()
    pd.DataFrame({'ID': ['ref1', 'ref2'], 'RHEA_ID': [10001, 10002]}).to_csv(
        tmp_path / 'rhea/rhea2uniprot_sprot.tsv', sep='\t', index=False)
    if kind == 'no_columns':
        producer = pd.DataFrame()
    elif kind == 'header_only':
        producer = pd.DataFrame(columns=[0, 1])
    else:
        producer = pd.DataFrame([['ref1', 'target1'], ['ref2', 'target2']])
    producer.to_excel(work / 'test_rhea_homolog0.9.xlsx')
    process_rhea('test', 'ecoli', .9, [])
    result = pd.read_excel(work / 'test_rheaid_plmsearch.xlsx', index_col=0)
    assert list(result.columns) == ['ID', 'Query', 'Rhea_id']
    assert result.values.tolist() == ([['ref1', 'target1', 10001], ['ref2', 'target2', 10002]]
                                      if kind == 'hits' else [])
