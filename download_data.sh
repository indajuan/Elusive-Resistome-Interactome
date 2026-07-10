#!/usr/bin/env bash
#
# Downloads dataset files for the Elusive_app project from Zenodo (record
# 19702877) and the IndaDiaz2026__ARGTools GitHub repository.
#
# Usage:
#   ./download_data.sh [destination_dir]
#
# destination_dir defaults to ./data

set -euo pipefail

DEST_DIR="${1:-data}"
mkdir -p "$DEST_DIR"

ZENODO_FILES=(
  "abundance_richness.csv.gz"
  "args_abundances.tsv.gz"
  "conversion_aro_geneclass.csv"
  "metagenomes_metadata.csv"
  "pipeline_abricate.argannot.norm.csv"
  "pipeline_abricate.card.norm.csv"
  "pipeline_abricate.megares.norm.csv"
  "pipeline_abricate.ncbi.norm.csv"
  "pipeline_abricate.resfinder.norm.csv"
  "pipeline_amrfinder.norm.csv"
  "pipeline_amrfinder.norm.prot.csv"
  "pipeline_deeparg.norm.csv"
  "pipeline_deeparg.norm.prot.csv"
  "pipeline_fargene.csv"
  "pipeline_fargene.prot.csv"
  "pipeline_resfinder.norm.csv"
  "pipeline_rgi.blast.csv"
  "pipeline_rgi.diamond.csv"
  "pipeline_rgi.diamond.prot.csv"
  "reported_unigenes_as_ARG_per_habitat.csv"
)

ZENODO_BASE="https://zenodo.org/records/19702877/files"
GITHUB_RAW_URL="https://raw.githubusercontent.com/BigDataBiology/IndaDiaz2026__ARGTools/main/cluster_vsearch/clusters.uc"

download() {
  local url="$1"
  local out="$2"
  if [[ -f "$out" ]]; then
    echo "Skipping (already exists): $out"
    return
  fi
  echo "Downloading: $out"
  curl -fL --retry 3 --retry-delay 5 -o "$out.part" "$url"
  mv "$out.part" "$out"
}

for fname in "${ZENODO_FILES[@]}"; do
  download "${ZENODO_BASE}/${fname}?download=1" "${DEST_DIR}/${fname}"
done

download "$GITHUB_RAW_URL" "${DEST_DIR}/clusters.uc"

echo "All downloads complete. Files saved in: $DEST_DIR"
