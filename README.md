# Elusive_app

Data pipeline and interactive web app companion to *"The elusive resistome: a
global comparison reveals large discrepancies among detection pipelines"*
(Inda-Díaz et al., bioRxiv 2026). Rebuilt from the paper's
Zenodo data record and VSEARCH clustering output (GitHub).

## Overview

```
download_data.sh    ->  data_zenodo_github/         (raw source files)
build_unigenes.py    ->  unigenes.tsv                (per-tool ARG calls, habitat-filtered)
rarefy_abundances.py ->  (rarefaction helper, used by core_resistome.py)
build_app_data.py    ->  webapp/data/*.json          (everything the web app needs, except core/pan)
core_resistome.py    ->  core_resistome.tsv           (core-resistome subsampling analysis)
serve_webapp.sh       ->  http://localhost:8010       (serves webapp/)
```

Run the steps in this order the first time; each script is safe to re-run
(downloads are skipped if the file already exists; the rest just overwrite
their output).

## 1. Download the raw data

```bash
./download_data.sh data_zenodo_github
```

Downloads 20 files from the Zenodo record (record 19702877) plus
`clusters.uc` (vsearch cluster assignments) from the
`BigDataBiology/IndaDiaz2026__ARGTools` GitHub repo, into
`data_zenodo_github/`. Includes the per-pipeline ARG call tables
(`pipeline_*.csv`), sample metadata (`metagenomes_metadata.csv`), gene-class
abundance/richness (`abundance_richness.csv.gz`), per-gene abundance
(`args_abundances.tsv.gz`), and the habitat-occurrence table
(`reported_unigenes_as_ARG_per_habitat.csv`).

```bash
./download_data.sh [destination_dir]   # destination_dir defaults to ./data
```

## 2. Build the unigenes table

```bash
./build_unigenes.py [data_dir] [output_tsv]
# data_dir defaults to data_zenodo_github, output_tsv defaults to unigenes.tsv
```

For each of the tool/identity-threshold variants (DeepARG, fARGene, the
five ABRicate databases, RGI, AMRFinder-Plus, ResFinder, plus their
70/80/90%-identity and amino-acid/BLAST variants), tags every ARG call with
its `tool` name, maps its vsearch cluster centroid, and derives
`new_level_centroid` (the ARG class, propagated across every tool sharing a
centroid via a global majority vote — see the script's diagnostics output
for clusters with disagreeing classifications). Rows are kept only if the
gene was reported as an ARG in at least one of 13 target habitats (human
gut/oral/skin/nose/vagina, dog/cat/mouse/pig gut, wastewater, marine,
freshwater, soil). Writes `unigenes.tsv` (~565k rows).

## 3. Build the web app data

```bash
./build_app_data.py [data_dir] [unigenes_tsv] [out_dir]
# defaults: data_zenodo_github, unigenes.tsv, webapp/data
```

Reads `unigenes.tsv` plus `abundance_richness.csv.gz` and
`metagenomes_metadata.csv`, and reproduces every JSON file the front-end
needs (ARG counts, Jaccard index, identity distributions, gene-class
breakdowns, class-specific coverage, habitat-resolved abundance/richness,
plus per-habitat variants of all of the above) — **except**
`pan_resistome.json`/`core_resistome.json`, which come from step 4 instead.
Also writes `table_s3_full.csv.gz` (the full per-gene ARG table, offered as
a direct download in the app).

## 4. Core-resistome analysis (optional, not yet wired into the app)

## 5. Run the web app

```bash
./serve_webapp.sh
```

Serves `webapp/` at **http://localhost:8010** (needs `webapp/data/*.json`
from step 3 to already exist). Open that URL in a browser — it's a static
site (HTML + JS + JSON, via Plotly), so any static file server works
equally well if you'd rather not use the provided script.

The app has a collapsible left-hand menu (Introduction / General analysis /
Habitat level, each with its own submenu) instead of a single scrolling
page — click any item to jump straight to that view.

## Requirements

- `curl` (step 1)
- Python 3 with `pandas` and `numpy` (steps 2–4)
- Any modern browser (step 5) — Plotly is loaded from a CDN, so an internet
  connection is needed even though the app itself runs entirely client-side

## To build

```bash
./core_resistome.py [data_dir] [output_tsv] [sub_sampling_size] [number_of_subsamples] [depth] [seed]
# defaults: data_zenodo_github, core_resistome.tsv, 100, 500, 5e6, 2000
```

Repeatedly (`number_of_subsamples` times) draws `sub_sampling_size` samples
per habitat without replacement, rarefies `args_abundances.tsv.gz` to
`depth` reads per sample (via `rarefy_abundances.py`, seeded from `seed`),
and tallies how consistently each ARG class qualifies as "core" (present in
≥p of a subsample's samples) across a range of thresholds. This is a
faithful-but-adapted port of the paper's R `core_resistome()`/
`filter_samples_core()` pipeline — see the script's docstring for the exact
differences from the original R code.

This is computationally the heaviest step (500 iterations × 21 tools over
the full abundance table) and hasn't been run end-to-end yet; the app's
**Pan- and Core-resistome** page currently shows a "coming soon" placeholder
pending that run and the matching `build_app_data.py`-style JSON export.

`rarefy_abundances.py` can also be run standalone to inspect the rarefied
abundance table on its own:

```bash
./rarefy_abundances.py [data_dir] [output_file] [depth] [seed]
# defaults: data_zenodo_github, data_zenodo_github/args_abundances_rarefied.tsv.gz, 5e6, 2000
```
