#!/usr/bin/env bash
# =============================================================================
# Downloads the Bulk2Spot tutorial dataset: NanoString's public WTA kidney demo
# (Diabetic Kidney Disease vs normal), the dataset used in Bioconductor's
# "Analyzing GeoMx-NGS RNA Expression Data with GeomxTools" vignette:
# 239 DCC files, the Hs WTA v1.0 PKC and the annotation workbook (~75 MB).
#
# Source: github.com/Nanostring-Biostats/GeoMxWorkflows, inst/extdata/WTA_NGS_Example/
#
# Usage: tests/tutorial/download_tutorial_data.sh [PARALLEL_DOWNLOADS]   (default 8)
# Re-running is safe: files already downloaded are skipped.
# =============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
DCC_DIR="${DATA_DIR}/dccs"
PKC_DIR="${DATA_DIR}/pkcs"
ANNO_DIR="${DATA_DIR}/annotation"
PARALLEL="${1:-8}"
RAW_BASE="https://raw.githubusercontent.com/Nanostring-Biostats/GeoMxWorkflows/main/inst/extdata/WTA_NGS_Example"
ANNO_FILE="kidney_AOI_Annotations_all_vignette.xlsx"

mkdir -p "${DCC_DIR}" "${PKC_DIR}" "${ANNO_DIR}"

# Download to a .tmp file first, so an interrupted download is never mistaken for a complete one
fetch() {
  local url="$1" out="$2"
  [[ -s "${out}" ]] && return 0
  curl -sSL --retry 3 "${url}" -o "${out}.tmp" && mv "${out}.tmp" "${out}"
}
export -f fetch

echo "== Annotation workbook =="
fetch "${RAW_BASE}/annotation/${ANNO_FILE}" "${ANNO_DIR}/${ANNO_FILE}"

echo "== PKC (Hs WTA v1.0) =="
if [[ ! -s "${PKC_DIR}/Hsa_WTA_v1.0.pkc" ]]; then
  fetch "${RAW_BASE}/pkcs/Hsa_WTA_v1.0.pkc.zip" "${PKC_DIR}/Hsa_WTA_v1.0.pkc.zip"
  unzip -o -q "${PKC_DIR}/Hsa_WTA_v1.0.pkc.zip" -d "${PKC_DIR}"
fi

echo "== DCC files (239 segments, ${PARALLEL} parallel downloads) =="
xargs -a "${SCRIPT_DIR}/dcc_filenames.txt" -P "${PARALLEL}" -I{} \
  bash -c 'fetch "$0/dccs/$1" "$2/$1"' "${RAW_BASE}" {} "${DCC_DIR}"

n_ok=$(find "${DCC_DIR}" -name "*.dcc" | wc -l)
echo "  ${n_ok}/239 DCC files present"
if (( n_ok < 239 )); then
  echo "  WARNING: some downloads failed -- re-run this script to retry." >&2
  exit 1
fi

echo
echo "Done. Raw data in: ${DATA_DIR}"
echo "Next (from the repository root): ./run.py -w report -c tests/tutorial/config_tutorial.yaml -q 8"
