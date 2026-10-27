#!/usr/bin/env bash
set -euo pipefail
BASE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"
PYTHON="${PYTHON:-python}"
if [ "${DEBUG:-0}" = "1" ]; then set -x; fi
log() { printf '[setup] %s\n' "$*"; }
trap 'echo "[setup] FAILED at line $LINENO"; exit 1' ERR
log "Starting setup in $BASE_DIR"
if [ -r /etc/os-release ]; then
  . /etc/os-release
  if [ "${ID:-}" != "ubuntu" ]; then
    echo "Ubuntu only. Detected: ${ID:-unknown}" >&2
    exit 1
  fi
  log "OS check passed: ${ID:-unknown} ${VERSION_ID:-}"
else
  echo "Ubuntu only. /etc/os-release not found." >&2
  exit 1
fi
if [ ! -f "resource-manifest.json" ] || [ "${ALPHAGEM_EXTRACT_RESOURCES:-0}" = "1" ]; then
  if [ -f "data_source.tar.xz" ]; then
    sha256sum -c data_source.tar.xz.sha256
    log "Extracting data_source.tar.xz"
    tar -xJf data_source.tar.xz
    log "Extraction complete"
  else
    echo "Put the matching data_source.tar.xz in this checkout first." >&2
    exit 1
  fi
else
  log "Resource manifest already present; verifying installed files"
fi
"$PYTHON" scripts/verify_resources.py
mkdir -p temp working struct_data/taryeast tools/CLEAN/app/data/inputs tools/CLEAN/app/data/esm_data tools/CLEAN/app/results/inputs
for reference in yeast ecoli human strco synechocystis; do
  if [ -f "data_available/${reference}.fasta" ]; then
    makeblastdb -in "data_available/${reference}.fasta" -dbtype prot
  fi
done
if ! command -v foldseek >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    log "Installing foldseek via conda"
    conda install -y -c conda-forge -c bioconda foldseek=9.427df8a
    log "foldseek install finished"
  else
    log "conda not found, skip foldseek install"
  fi
else
  log "foldseek already present"
fi
if ! command -v diamond >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    log "Installing diamond via conda"
    conda install -y -c conda-forge -c bioconda diamond=2.1.10
    log "diamond install finished"
  else
    log "conda not found, skip diamond install"
  fi
else
  log "diamond already present"
fi
if [ "${ALPHAGEM_SKIP_TOOL_DATABASES:-0}" = "1" ]; then
  log "Tool databases deferred; install them before running annotation"
else
  read -r -a eggnog_taxa <<< "${ALPHAGEM_EGGNOG_TAXA:-fungi}"
  "$PYTHON" scripts/install_tool_databases.py --taxa "${eggnog_taxa[@]}"
fi
if [ -d "./tools/CLEAN/app" ]; then
  log "Installing CLEAN from ./tools/CLEAN/app"
  pushd ./tools/CLEAN/app >/dev/null
  "$PYTHON" build.py install
  popd >/dev/null
  if "$PYTHON" - <<'PY'
import sys, importlib.util
sys.exit(0 if (importlib.util.find_spec('CLEAN') or importlib.util.find_spec('clean')) else 1)
PY
  then
    log "CLEAN import available"
  else
    log "CLEAN import not found after install"
  fi
else
  log "CLEAN directory not found, skip"
fi
if [ -d "./tools/OrthoFinder" ]; then
  log "Preparing OrthoFinder"
  pushd ./tools/OrthoFinder >/dev/null
  "$PYTHON" setup.py install
  popd >/dev/null
else
  log "OrthoFinder directory not found, skip"
fi
log "Setup completed successfully"
