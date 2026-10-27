"""Gapfill genes must use complete tokens and only fresh, successful alignments."""

from pathlib import Path
import subprocess
import sys
from unittest.mock import Mock

import cobra
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src')]
from src import gapfilling


@pytest.fixture
def work(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / 'working/test'
    work.mkdir(parents=True)
    (tmp_path / 'data_available').mkdir()
    pd.DataFrame([
        ['ref1', 'ref1', 'YER060W'], ['ref2', 'ref2', 'YER060W-A'],
        ['ref3', 'ref3', 'OTHER'],
    ], columns=['Entry', 'Entry Name', 'Gene Names (ordered locus)']).to_excel(
        tmp_path / 'data_available/yeast.xlsx', index=False)
    return work


def hit(reference, target, score):
    return '\t'.join([reference, target, '100', '100', '0', '0', '1', '100', '1', '100', '0', str(score)]) + '\n'


def model():
    result = cobra.Model('regression')
    result.solver = 'glpk'
    reaction = cobra.Reaction('R1')
    reaction.gene_reaction_rule = 'YER060W or (YER060W-A and OTHER)'
    result.add_reactions([reaction])
    return result


def test_full_token_replacement_and_highest_score(work):
    (work / 'outfiletest.csv').write_text(
        '# BLASTP result\n' + hit('sp|ref1|ref1', 'sp|weak|weak', 200)
        + hit('sp|ref1|ref1', 'sp|target|target', 300)
        + hit('sp|ref3|ref3', 'sp|low|low', 199))
    result, unmatched = gapfilling.read_gapgenes(
        model(), ['sp|ref1|ref1', 'sp|ref3|ref3'], ['R1'], 'yeast', 'test')
    assert result.reactions.R1.gpr.genes == {'target', 'YER060W-A', 'OTHER'}
    assert {g.id for g in result.genes} == {'target', 'YER060W-A', 'OTHER'}
    assert unmatched == ['OTHER']


@pytest.mark.parametrize('content', ['', '# BLASTP\n# 0 hits found\n'])
def test_empty_search_is_valid_absent_evidence(work, content):
    (work / 'outfiletest.csv').write_text(content)
    original = model()
    _, unmatched = gapfilling.read_gapgenes(original, ['sp|ref1|ref1'], ['R1'], 'yeast', 'test')
    assert unmatched == ['YER060W']
    assert original.reactions.R1.gpr.genes == {'YER060W', 'YER060W-A', 'OTHER'}


@pytest.mark.parametrize('content', ['truncated\trow\n', hit('ref1', 'target', 'bad'),
                                    hit('ref1', 'target', 'inf')])
def test_malformed_results_fail_instead_of_silently_dropping_genes(work, content):
    (work / 'outfiletest.csv').write_text(content)
    with pytest.raises(ValueError):
        gapfilling.read_gapgenes(model(), ['sp|ref1|ref1'], ['R1'], 'yeast', 'test')


def test_missing_mapping_is_explicit(work):
    (work / 'outfiletest.csv').write_text('')
    with pytest.raises(ValueError, match='missing ordered locus mapping: absent'):
        gapfilling.read_gapgenes(model(), ['sp|absent|absent'], ['R1'], 'yeast', 'test')


@pytest.mark.parametrize('nucleotide', [False, True])
@pytest.mark.parametrize('fail_at', [1, 2])
def test_external_failure_stops_without_publishing_partial_output(work, monkeypatch, nucleotide, fail_at):
    if nucleotide:
        (work / 'test_non_anno.fasta').touch()
    output = work / 'outfiletest.csv'
    output.write_text('previous alignment')
    calls = []

    def run(command, check):
        assert check is True
        calls.append(command)
        if len(calls) == fail_at:
            raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(gapfilling.subprocess, 'run', run)
    with pytest.raises(subprocess.CalledProcessError):
        gapfilling.run_gapgene_alignment('test', 'yeast')
    assert len(calls) == fail_at
    assert output.read_text() == 'previous alignment'
    assert not list(work.glob('gapgene-*'))


@pytest.mark.parametrize('nucleotide', [False, True])
def test_success_publishes_fresh_output(work, monkeypatch, nucleotide):
    if nucleotide:
        (work / 'test_non_anno.fasta').touch()
    calls = []

    def run(command, check):
        calls.append(command)
        if len(calls) == 2:
            option = '-out' if nucleotide else '--out'
            Path(command[command.index(option) + 1]).write_text('# fresh, no hits\n')

    monkeypatch.setattr(gapfilling.subprocess, 'run', run)
    gapfilling.run_gapgene_alignment('test', 'yeast')
    assert calls[0][0] == ('makeblastdb' if nucleotide else 'diamond')
    assert calls[1][0] == ('tblastn' if nucleotide else 'diamond')
    assert (work / 'outfiletest.csv').read_text() == '# fresh, no hits\n'


def test_success_without_output_is_an_error(work, monkeypatch):
    monkeypatch.setattr(gapfilling.subprocess, 'run', Mock())
    with pytest.raises(RuntimeError, match='without producing'):
        gapfilling.run_gapgene_alignment('test', 'yeast')
