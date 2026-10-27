"""Cache provenance must follow sequences, model weights and serialized output."""

import json
from pathlib import Path
import sys
from unittest.mock import Mock

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import plmsearch_embedding as embedding


@pytest.fixture
def cache(tmp_path, monkeypatch):
    fasta = tmp_path / 'target.fasta'
    fasta.write_text('>sp|same_id|same_id\nMAAA\n')
    checkpoint = tmp_path / 'model.pt'
    checkpoint.write_bytes(b'weights-v1')
    output = tmp_path / 'target_embedding.pkl'
    compute = Mock(return_value={'sp|same_id|same_id': torch.ones(3)})
    monkeypatch.setattr(embedding, '_compute_embeddings', compute)
    return fasta, checkpoint, output, compute


def generate(cache):
    fasta, checkpoint, output, _ = cache
    return embedding._generate_cached_embeddings(fasta, output, checkpoint, True)


def test_same_inputs_reuse_without_inference(cache):
    assert generate(cache)['status'] == 'generated'
    assert generate(cache)['status'] == 'reused'
    cache[3].assert_called_once()


@pytest.mark.parametrize('changed', ['sequence', 'checkpoint', 'recipe', 'output',
                                    'missing_manifest', 'broken_manifest', 'legacy_manifest',
                                    'incomplete_manifest', 'missing_output'])
def test_changed_inputs_or_unverified_output_regenerate(cache, monkeypatch, changed):
    fasta, checkpoint, output, compute = cache
    generate(cache)
    manifest = output.with_suffix('.manifest.json')
    if changed == 'sequence':
        fasta.write_text('>sp|same_id|same_id\nMCCC\n')
    elif changed == 'checkpoint':
        checkpoint.write_bytes(b'weights-v2')
    elif changed == 'recipe':
        monkeypatch.setitem(embedding.EMBEDDING_RECIPE, 'version', 'next-recipe')
    elif changed == 'output':
        output.write_bytes(b'corrupt pickle')
    elif changed == 'missing_manifest':
        manifest.unlink()
    elif changed == 'broken_manifest':
        manifest.write_text('{broken json')
    elif changed == 'legacy_manifest':
        manifest.write_text('[]')
    elif changed == 'incomplete_manifest':
        data = json.loads(manifest.read_text())
        data['status'] = 'running'
        manifest.write_text(json.dumps(data))
    else:
        output.unlink()
    assert generate(cache)['status'] == 'generated'
    assert compute.call_count == 2
    assert generate(cache)['status'] == 'reused'


def test_generation_failure_preserves_previous_pair(cache):
    generate(cache)
    fasta, _, output, compute = cache
    manifest = output.with_suffix('.manifest.json')
    previous = output.read_bytes(), manifest.read_bytes()
    fasta.write_text('>sp|same_id|same_id\nMCCC\n')
    compute.side_effect = RuntimeError('inference failed')
    with pytest.raises(RuntimeError, match='inference failed'):
        generate(cache)
    assert previous == (output.read_bytes(), manifest.read_bytes())
    assert not embedding._cache_matches(output, manifest, embedding._input_fingerprint(fasta, cache[1]))


@pytest.mark.parametrize('changed', ['sequence', 'checkpoint'])
def test_input_mutation_during_generation_is_not_published(cache, changed):
    fasta, checkpoint, output, compute = cache

    def mutate(*args):
        (fasta if changed == 'sequence' else checkpoint).write_bytes(b'changed')
        return {'id': torch.ones(3)}

    compute.side_effect = mutate
    with pytest.raises(RuntimeError, match='changed during'):
        generate(cache)
    assert not output.exists()
    assert not output.with_suffix('.manifest.json').exists()


def test_named_model_first_download_is_recorded(tmp_path, monkeypatch):
    monkeypatch.setattr(torch.hub, 'get_dir', lambda: str(tmp_path / 'hub'))
    fasta = tmp_path / 'test.fasta'
    fasta.write_text('>id\nMAAA\n')
    checkpoint = tmp_path / 'hub/checkpoints/esm1b_t33_650M_UR50S.pt'

    def download(*args):
        checkpoint.parent.mkdir(parents=True)
        checkpoint.write_bytes(b'downloaded weights')
        return {'id': torch.ones(3)}

    compute = Mock(side_effect=download)
    monkeypatch.setattr(embedding, '_compute_embeddings', compute)
    output = tmp_path / 'test_embedding.pkl'
    assert embedding._generate_cached_embeddings(fasta, output, 'esm1b_t33_650M_UR50S', True)['status'] == 'generated'
    assert embedding._generate_cached_embeddings(fasta, output, 'esm1b_t33_650M_UR50S', True)['status'] == 'reused'
    manifest = json.loads(output.with_suffix('.manifest.json').read_text())
    assert manifest['inputs']['model']['checkpoint_sha256'] == embedding._sha256(checkpoint)
    compute.assert_called_once()


def test_empty_fasta_fails_before_loading_model(tmp_path, monkeypatch):
    fasta = tmp_path / 'empty.fasta'
    fasta.touch()
    load = Mock()
    monkeypatch.setattr(embedding.pretrained, 'load_model_and_alphabet', load)
    with pytest.raises(ValueError, match='empty FASTA'):
        embedding._compute_embeddings(fasta, 'model', True)
    load.assert_not_called()
