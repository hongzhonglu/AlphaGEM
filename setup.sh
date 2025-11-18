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
if [ ! -d "./tools" ] || [ ! -d "./data_available" ] || [ ! -d "./plmsearchtools" ] || [ ! -d "./struct_data" ]; then
  if [ -f "data_source.tar.xz" ]; then
    log "Extracting data_source.tar.xz"
    tar -xJf data_source.tar.xz
    log "Extraction complete"
  else
    log "data_source.tar.xz not found, skip extraction"
  fi
else
  log "Data directories already present, skip extraction"
fi
if ! command -v foldseek >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    log "Installing foldseek via conda"
    conda install -y -c bioconda -c conda-forge foldseek
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
    conda install -y -c bioconda -c conda-forge diamond
    log "diamond install finished"
  else
    log "conda not found, skip diamond install"
  fi
else
  log "diamond already present"
fi
if [ -d "./tools/eggnog_mapper" ]; then
  if [ ! -d "./tools/eggnog_mapper/data" ] || [ ! -f "./tools/eggnog_mapper/data/eggnog.db" ]; then
    log "Downloading eggNOG data"
    pushd ./tools/eggnog_mapper >/dev/null
    yes | "$PYTHON" download_eggnog_data.py
    popd >/dev/null
    log "eggNOG data ready"
  else
    log "eggNOG data already present"
  fi
else
  log "eggnog_mapper directory not found, skip"
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
