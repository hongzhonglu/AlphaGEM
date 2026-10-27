# Runtime resources

Use the resource bundle matching `data_source.tar.xz.sha256` in this checkout.
The previous Figshare bundle has a different checksum; the updated download link
will be added after the new file is uploaded. Do not substitute the old bundle
for this snapshot.

```bash
# From the AlphaGEM repository root, with the downloaded archive here:
sha256sum -c data_source.tar.xz.sha256
tar -xJf data_source.tar.xz
conda run -n AlphaGEM python scripts/verify_resources.py
```

The archive contains the original Rhea resources, PLM model, matching reference
FASTA/locus tables/embeddings, inference code and trained checkpoints, and the
four reference structure snapshots (`ecoli`, `human`, `strco`, `yeast`). It includes
the corrected E. coli accession mappings to `b3036` and `b4223`.

The existing Git-tracked `models/reaction_pool_final.json` and the formal legacy
`general2specific.pkl` snapshot are retained. No new reaction pool, SQLite
annotation catalogue, or specificity backend is part of this update.

Target-species PDBs, historical reconstruction outputs, test/example datasets,
build caches and generated BLAST/Foldseek indexes are omitted. Supply your own
target FASTA and structures. `synechocystis` includes sequence/locus/embedding
resources for PLM mode; no matching structure snapshot is supplied. `universal`
still lacks a locus table and is not a validated complete workflow.

## Downloadable tool databases

eggNOG-mapper 2.1.12 / database 5.0.2 databases and its full source FASTA are
deliberately absent. Use its own download and taxon-selection/build methods:

```bash
# Preview the methods without downloading anything:
conda run -n AlphaGEM python scripts/install_tool_databases.py --dry-run --taxa fungi

# yeast reference:
conda run --no-capture-output -n AlphaGEM python scripts/install_tool_databases.py --taxa fungi
# ecoli / strco / synechocystis references:
conda run --no-capture-output -n AlphaGEM python scripts/install_tool_databases.py --taxa bacteria
# human reference:
conda run --no-capture-output -n AlphaGEM python scripts/install_tool_databases.py --taxa mammalia
# Or prepare all three groups in one run:
conda run --no-capture-output -n AlphaGEM python scripts/install_tool_databases.py --taxa fungi bacteria mammalia
```

The installer calls the bundled tool's annotation/taxonomy downloads, full
protein/taxon downloads, taxon selection, FASTA filtering and DIAMOND builder.
It uses HTTPS, propagates command failures, checks readable SQLite schemas and
native DIAMOND databases, and publishes a taxon database only after a successful
build. Downloading and building these databases needs substantial network and
disk capacity; the main annotation database alone is about 41 GB uncompressed.
New downloads were not performed as part of the offline packaging validation.

The complete eggNOG source FASTA is pinned by SHA-256. If an old extraction
contains the tiny partial `e5.proteomes.faa`, the installer stops explicitly;
move that file aside and rerun so the tool downloads the complete snapshot.
Do not treat the existence of a FASTA or `.dmnd` file as proof that installation
completed.

The unused Hugging Face `prot_bert_bfd/pytorch_model.bin` is omitted. DeepProZyme
uses its dedicated trained `model/model.pth`; matching tokenizer files remain.
Its small fixed `swissprot_enzyme_diamond.dmnd` reference is retained because it
is a tool-specific inference input, with no equivalent pinned downloader in the
bundled tool. CLEAN weights, calibration data and the project-curated Rhea
embeddings are also retained.

## Setup and indexes

```bash
# Verify/extract resources, rebuild reference BLAST indexes, install tools,
# and install the fungal eggNOG databases by default:
conda run --no-capture-output -n AlphaGEM bash setup.sh

# Prepare resources/tools now, defer large database downloads:
ALPHAGEM_SKIP_TOOL_DATABASES=1 conda run --no-capture-output -n AlphaGEM bash setup.sh

# Choose different eggNOG groups during setup:
ALPHAGEM_EGGNOG_TAXA='bacteria mammalia' conda run --no-capture-output -n AlphaGEM bash setup.sh
```

`ALPHAGEM_EXTRACT_RESOURCES=1` forces extraction of the matching archive when
upgrading an existing installation. Extraction does not delete existing tool
databases or user outputs. The embedded manifest verifies every supplied file;
downloaded/generated files are outside that manifest. Target ESM weights and
ProtT5 weights are obtained by their normal model-loading methods when needed.
