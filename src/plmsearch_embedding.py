import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import torch
import pickle
from esm.data import FastaBatchedDataset
from esm import pretrained
import torch.serialization as ts

if hasattr(ts, "add_safe_globals"):
    ts.add_safe_globals([argparse.Namespace])

EMBEDDING_RECIPE = {
    "version": "esm1b_legacy_mean_layer33_v1",
    "representation_layer": 33,
    "truncation_seq_length": 1022,
    "token_budget": 16384,
}


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _input_fingerprint(fasta, esm_model_path):
    model_name = str(esm_model_path)
    checkpoint = (Path(model_name) if model_name.endswith(".pt") else
                  Path(torch.hub.get_dir()) / "checkpoints" / f"{model_name}.pt")
    return {
        "fasta_sha256": _sha256(fasta),
        "model": {
            "identifier": model_name,
            "checkpoint_sha256": _sha256(checkpoint) if checkpoint.is_file() else None,
        },
        "recipe": dict(EMBEDDING_RECIPE),
    }


def _cache_matches(result, manifest_path, fingerprint):
    if not result.is_file() or not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return (isinstance(manifest, dict) and manifest.get("schema_version") == 1
                and manifest.get("status") == "completed"
                and manifest.get("inputs") == fingerprint
                and manifest.get("embedding_sha256") == _sha256(result))
    except (OSError, ValueError):
        return False


def _compute_embeddings(fasta, esm_model_path, nogpu):
    dataset = FastaBatchedDataset.from_file(str(fasta))
    if not len(dataset):
        raise ValueError(f"Cannot generate embeddings from an empty FASTA: {fasta}")
    esm_model, alphabet = pretrained.load_model_and_alphabet(str(esm_model_path))
    esm_model.eval()
    use_gpu = torch.cuda.is_available() and not nogpu
    if use_gpu:
        esm_model = esm_model.cuda()
        print("Transferred model to GPU")
    batches = dataset.get_batch_indices(EMBEDDING_RECIPE["token_budget"], extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset,
        collate_fn=alphabet.get_batch_converter(EMBEDDING_RECIPE["truncation_seq_length"]),
        batch_sampler=batches,
    )
    print(f"Read {fasta} with {len(dataset)} sequences")
    embeddings = {}
    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
            if use_gpu:
                toks = toks.to(device="cuda", non_blocking=True)
            out = esm_model(toks, repr_layers=[33], return_contacts=False)["representations"][33]
            for index, label in enumerate(labels):
                # Preserve the existing pooling/truncation behavior.
                embeddings[label] = out[index, 1:len(strs[index]) + 1].mean(0).clone().cpu()
    return embeddings


def _generate_cached_embeddings(fasta, result, esm_model_path, nogpu):
    fasta, result = Path(fasta), Path(result)
    manifest_path = result.with_suffix(".manifest.json")
    fingerprint = _input_fingerprint(fasta, esm_model_path)
    if _cache_matches(result, manifest_path, fingerprint):
        print(f"Reusing verified ESM embeddings: {result}")
        return {"status": "reused", "manifest_path": str(manifest_path)}
    print(f"Generating ESM embeddings; cache absent, unverified or mismatched: {result}")
    embeddings = _compute_embeddings(fasta, esm_model_path, nogpu)
    current = _input_fingerprint(fasta, esm_model_path)
    if (current["fasta_sha256"] != fingerprint["fasta_sha256"] or
            (fingerprint["model"]["checkpoint_sha256"] is not None
             and current["model"] != fingerprint["model"])):
        raise RuntimeError("FASTA or ESM checkpoint changed during embedding generation; no cache written")
    result.parent.mkdir(parents=True, exist_ok=True)
    temporary_paths = []
    try:
        with tempfile.NamedTemporaryFile(dir=result.parent, suffix=".pkl.tmp", delete=False) as handle:
            temporary_result = Path(handle.name)
            temporary_paths.append(temporary_result)
            pickle.dump(embeddings, handle, protocol=pickle.HIGHEST_PROTOCOL)
        manifest = {
            "schema_version": 1, "status": "completed", "inputs": current,
            "fasta_path": str(fasta), "sequence_count": len(embeddings),
            "embedding_sha256": _sha256(temporary_result),
        }
        with tempfile.NamedTemporaryFile(
                dir=result.parent, mode="w", encoding="utf-8", suffix=".json.tmp", delete=False) as handle:
            temporary_manifest = Path(handle.name)
            temporary_paths.append(temporary_manifest)
            json.dump(manifest, handle, indent=2)
            handle.write("\n")
        temporary_result.replace(result)
        temporary_manifest.replace(manifest_path)
    finally:
        for path in temporary_paths:
            path.unlink(missing_ok=True)
    return {"status": "generated", "manifest_path": str(manifest_path)}


def embedding_generate(name, esm_model_path="esm1b_t33_650M_UR50S", nogpu=False):
    work = Path("working") / name
    return _generate_cached_embeddings(
        work / f"{name}.fasta", work / f"{name}_embedding.pkl", esm_model_path, nogpu,
    )


def embedding_generate_base(name, esm_model_path="esm1b_t33_650M_UR50S", nogpu=False):
    return _generate_cached_embeddings(
        Path("data_available") / f"{name}.fasta",
        Path("working") / f"{name}_embedding.pkl", esm_model_path, nogpu,
    )
