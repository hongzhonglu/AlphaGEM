"""Prevent substituting a differently normalized structural similarity score."""

from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import foldseek_knock1


def test_query_normalized_score_is_not_a_fallback(monkeypatch):
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        text = ('query,target,qtmscore,prob,qcov,tcov' if '-h' in command
                else '1.3c64211\n')
        return subprocess.CompletedProcess(command, 0, text, '')

    monkeypatch.setattr(foldseek_knock1.subprocess, 'run', run)
    with pytest.raises(RuntimeError, match='alignment-normalized TM scores'):
        foldseek_knock1.foldseekfind('/unused', 'unused', 'unused')
    assert calls == [['foldseek', 'convertalis', '-h'], ['foldseek', 'version']]


def test_required_fields_accept_supported_binary(monkeypatch):
    def run(command, **kwargs):
        assert command == ['foldseek', 'convertalis', '-h']
        return subprocess.CompletedProcess(
            command, 0, 'query,target,alntmscore,qtmscore,prob,qcov,tcov', '')

    monkeypatch.setattr(foldseek_knock1.subprocess, 'run', run)
    foldseek_knock1._require_foldseek_output()
