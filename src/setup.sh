#!/usr/bin/env bash
set -euo pipefail
BASE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"
PYTHON="${PYTHON:-python}"
if [ -r /etc/os-release ]; then
  . /etc/os-release
  if [ "${ID:-}" != "ubuntu" ]; then
    echo "Ubuntu only. Detected: ${ID:-unknown}" >&2
    exit 1
  fi
else
  echo "Ubuntu only. /etc/os-release not found." >&2
  exit 1
fi
if [ ! -d "./tools" ] || [ ! -d "./data_available" ] || [ ! -d "./plmsearchtools" ] || [ ! -d "./struct_data" ]; then
  if [ -f "data_source.tar.xz" ]; then
    tar -xJf data_source.tar.xz
  fi
fi
if ! command -v foldseek >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    conda install -y -c bioconda -c conda-forge foldseek
  fi
fi
if ! command -v diamond >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    conda install -y -c bioconda -c conda-forge diamond
  fi
fi
if [ -d "./tools/eggnog_mapper" ]; then
  if [ ! -d "./tools/eggnog_mapper/data" ] || [ ! -f "./tools/eggnog_mapper/data/eggnog.db" ]; then
    pushd ./tools/eggnog_mapper >/dev/null
    yes | "$PYTHON" download_eggnog_data.py
    popd >/dev/null
  fi
fi
if [ -d "./tools/clean" ]; then
  pushd ./tools/clean >/dev/null
  "$PYTHON" build.py install
  popd >/dev/null
fi
if [ -d "./tools/OrthoFinder" ]; then
  pushd ./tools/OrthoFinder >/dev/null
  "$PYTHON" setup.py
  popd >/dev/null
fi
