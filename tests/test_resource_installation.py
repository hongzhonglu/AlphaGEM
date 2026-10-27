"""Check upload-package integrity and refuse partial eggNOG source databases."""

import hashlib
import json
from pathlib import Path
import sqlite3
import sys
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import install_tool_databases as installer
from scripts.verify_resources import verify


def test_manifest_verification_checks_content_not_just_size(tmp_path):
    resource = tmp_path / 'reference.fasta'
    resource.write_bytes(b'MAAA')
    manifest = {'schema_version': 1, 'file_count': 1, 'files': [
        {'path': resource.name, 'bytes': 4, 'sha256': hashlib.sha256(b'MAAA').hexdigest()},
    ]}
    (tmp_path / 'resource-manifest.json').write_text(json.dumps(manifest))
    assert verify(tmp_path) == 1
    resource.write_bytes(b'MCCC')
    with pytest.raises(ValueError, match='changed resource'):
        verify(tmp_path)


def test_manifest_cannot_verify_files_outside_checkout(tmp_path):
    (tmp_path / 'resource-manifest.json').write_text(json.dumps({
        'schema_version': 1, 'file_count': 1,
        'files': [{'path': '../outside', 'bytes': 0, 'sha256': ''}],
    }))
    with pytest.raises(ValueError, match='Unsafe'):
        verify(tmp_path)


def installed_tool(tmp_path, monkeypatch):
    tool = tmp_path / 'tools/eggnog_mapper'
    data = tool / 'data'
    data.mkdir(parents=True)
    (tool / 'download_eggnog_data.py').touch()
    for name in ['eggnog.db', 'eggnog.taxa.db']:
        with sqlite3.connect(data / name) as db:
            db.execute('CREATE TABLE reference (id TEXT)')
    download = SimpleNamespace(download_annotations=Mock(), download_taxa=Mock())
    create = SimpleNamespace(download_proteins=Mock(), download_taxid_info=Mock(),
                             parse_taxa_table=Mock(), parse_proteins=Mock(), create_diamond_db=Mock())
    monkeypatch.setattr(installer, 'load_tool', lambda tool, name, data:
                        download if name == 'download_eggnog_data' else create)
    monkeypatch.setitem(sys.modules, 'eggnogmapper.version',
                        SimpleNamespace(__VERSION__='2.1.12', __DB_VERSION__='5.0.2'))
    return data, create


def test_old_partial_fasta_is_never_used_for_database_build(tmp_path, monkeypatch):
    data, create = installed_tool(tmp_path, monkeypatch)
    (data / 'e5.proteomes.faa').write_text('>a partial old bundle\nMAAA\n')
    with pytest.raises(ValueError, match='not the pinned complete snapshot'):
        installer.install(tmp_path, ['fungi'])
    create.create_diamond_db.assert_not_called()
    create.parse_proteins.assert_not_called()
    assert not (data / 'fungi.dmnd').exists()


def test_existing_native_database_is_checked_before_reuse(tmp_path, monkeypatch):
    data, _ = installed_tool(tmp_path, monkeypatch)
    (data / 'fungi.dmnd').write_bytes(b'invalid database')
    check = Mock(side_effect=ValueError('invalid native database'))
    monkeypatch.setattr(installer, 'check_diamond', check)
    with pytest.raises(ValueError, match='invalid native database'):
        installer.install(tmp_path, ['fungi'])
    check.assert_called_once_with(data / 'fungi.dmnd')


def test_failed_build_does_not_publish_a_partial_taxon_database(tmp_path, monkeypatch):
    data, create = installed_tool(tmp_path, monkeypatch)
    (data / 'e5.proteomes.faa').write_bytes(b'valid fixture')
    (data / 'e5.taxid_info.tsv').write_bytes(b'fixture')
    monkeypatch.setattr(installer, 'sha256', lambda path: installer.PROTEINS_SHA256)
    create.parse_taxa_table.return_value = {'1'}
    create.parse_proteins.side_effect = lambda output, source, taxa: Path(output).write_text('>id\nMAAA\n')
    create.create_diamond_db.side_effect = RuntimeError('native builder failed')
    with pytest.raises(RuntimeError, match='native builder failed'):
        installer.install(tmp_path, ['fungi'])
    assert not (data / 'fungi.dmnd').exists()
    assert not list(data.glob('fungi-build-*'))
