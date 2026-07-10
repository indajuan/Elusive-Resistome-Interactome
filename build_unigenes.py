#!/usr/bin/env python3
"""
Builds unigenes.tsv by, for each ARG-detection pipeline output file:
  1. adding a "tool" column identifying the tool/run (and identity cutoff,
     where applicable),
  2. extracting the columns: query, tool, ARO, parent, parent_description,
     new_level, id,
  3. mapping each query to its vsearch cluster centroid (clusters.uc) and
     computing, per centroid, the centroid's own annotation and the
     majority-vote annotation across cluster members (diagnostics only —
     mirrors the R analysis, which computes these but does not fold them
     back into the kept columns),
then concatenating all pipelines together and keeping only rows whose
"query" (a unigene id) was reported as an ARG in one of a fixed set of
habitats, per reported_unigenes_as_ARG_per_habitat.csv.

Usage:
  ./build_unigenes.py [data_dir] [output_tsv]

data_dir defaults to data_zenodo_github, output_tsv defaults to unigenes.tsv
"""

import sys
from pathlib import Path

import pandas as pd

DATA_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data_zenodo_github")
OUTPUT_TSV = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("unigenes.tsv")

# (source file, tool label, minimum "id" (% identity) to keep, or None)
JOBS = [
    ("pipeline_deeparg.norm.csv", "DeepARG", None),
    ("pipeline_fargene.csv", "fARGene", None),
    ("pipeline_abricate.argannot.norm.csv", "ABRicate-ARGANNOT", None),
    ("pipeline_abricate.megares.norm.csv", "ABRicate-MEGARes", None),
    ("pipeline_rgi.diamond.csv", "RGI", None),
    ("pipeline_abricate.card.norm.csv", "ABRicate-CARD", None),
    ("pipeline_amrfinder.norm.prot.csv", "AMRFinder-Plus", None),
    ("pipeline_abricate.ncbi.norm.csv", "ABRicate-NCBI", None),
    ("pipeline_resfinder.norm.csv", "ResFinder", None),
    ("pipeline_abricate.resfinder.norm.csv", "ABRicate-ResFinder", None),
    ("pipeline_deeparg.norm.csv", "DeepARG-70%", 70),
    ("pipeline_deeparg.norm.csv", "DeepARG-80%", 80),
    ("pipeline_deeparg.norm.csv", "DeepARG-90%", 90),
    ("pipeline_rgi.diamond.csv", "RGI-70%", 70),
    ("pipeline_rgi.diamond.csv", "RGI-80%", 80),
    ("pipeline_rgi.diamond.csv", "RGI-90%", 90),
    ("pipeline_deeparg.norm.prot.csv", "DeepARG-aa", None),
    ("pipeline_rgi.blast.csv", "RGI-BLAST", None),
    ("pipeline_rgi.diamond.prot.csv", "RGI-aa", None),
    ("pipeline_fargene.prot.csv", "fARGene-aa", None),
    ("pipeline_amrfinder.norm.csv", "AMRFinder-Plus-nt", None),
]

OUTPUT_COLUMNS = ["query", "tool", "ARO", "parent", "parent_description", "new_level", "id"]

HABITATS = {
    "human gut", "human oral", "human skin", "human nose", "human vagina",
    "dog gut", "cat gut", "mouse gut", "pig gut",
    "wastewater", "marine", "freshwater", "soil",
}


def load_clusters():
    """clusters.uc (vsearch UC format) -> {query_id: centroid_id}."""
    cols = [f"V{i}" for i in range(1, 11)]
    df = pd.read_csv(DATA_DIR / "clusters.uc", sep="\t", header=None, names=cols, dtype=str)
    df = df[df["V1"] != "C"]  # keep only S(eed)/H(it) records
    df["centroid"] = df["V10"].where(df["V10"] != "*", df["V9"])
    return dict(zip(df["V9"], df["centroid"]))


def add_centroid_diagnostics(df: pd.DataFrame, centroid_map: dict, tool: str) -> pd.DataFrame:
    """Mirrors the R centroid / majority-rule analysis. Keeps new_level_centroid
    (the centroid's own annotation, falling back to the cluster majority vote)
    but drops the other helper columns (centroid, new_level_majority)."""
    df = df.copy()
    df["centroid"] = df["query"].map(centroid_map)

    n_distinct = df.groupby("centroid")["new_level"].transform(lambda s: s.nunique(dropna=False))
    heterogeneous = df.loc[n_distinct > 1, "centroid"].nunique()
    if heterogeneous:
        print(f"    {tool}: {heterogeneous:,} centroids have >1 distinct new_level among members")

    own_level = (
        df.drop_duplicates(subset="query", keep="first")
        .set_index("query")["new_level"]
        .to_dict()
    )
    df["new_level_centroid"] = df["centroid"].map(own_level)

    majority = (
        df.groupby("centroid")["new_level"]
        .agg(lambda s: s.dropna().mode().iloc[0] if not s.dropna().mode().empty else pd.NA)
        .to_dict()
    )
    df["new_level_majority"] = df["centroid"].map(majority)
    df["new_level_centroid"] = df["new_level_centroid"].fillna(df["new_level_majority"])

    return df.drop(columns=["centroid", "new_level_majority"])


def load_and_extract(fname: str, tool: str, min_id, centroid_map: dict):
    path = DATA_DIR / fname
    df = pd.read_csv(path, dtype=str, low_memory=False)

    # every pipeline file carries the resistance-class label under "geneclass_argcompare"
    df["new_level"] = df["geneclass_argcompare"]

    df["id"] = pd.to_numeric(df["id"], errors="coerce")
    if min_id is not None:
        df = df[df["id"] >= min_id]

    df["tool"] = tool

    for col in OUTPUT_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA

    extracted = df[OUTPUT_COLUMNS]
    return add_centroid_diagnostics(extracted, centroid_map, tool)


def allowed_queries():
    habitat_df = pd.read_csv(DATA_DIR / "reported_unigenes_as_ARG_per_habitat.csv", dtype=str)
    habitat_df = habitat_df[habitat_df["habitat"].isin(HABITATS)]
    return set(habitat_df["unigene"].unique())


def main():
    print("Loading clusters.uc ...")
    centroid_map = load_clusters()
    print(f"  {len(centroid_map):,} query -> centroid mappings")

    print(f"Reading habitat-filtered unigene set from reported_unigenes_as_ARG_per_habitat.csv ...")
    keep_queries = allowed_queries()
    print(f"  {len(keep_queries):,} unigenes qualify (reported in {sorted(HABITATS)})")

    parts = []
    for fname, tool, min_id in JOBS:
        label = f"{fname} -> tool={tool!r}" + (f", id>={min_id}" if min_id else "")
        print(f"Processing {label}")
        extracted = load_and_extract(fname, tool, min_id, centroid_map)
        before = len(extracted)
        extracted = extracted[extracted["query"].isin(keep_queries)]
        print(f"  {before:,} rows -> {len(extracted):,} after habitat filter")
        parts.append(extracted)

    merged = pd.concat(parts, ignore_index=True)
    merged.to_csv(OUTPUT_TSV, sep="\t", index=False)
    print(f"Wrote {len(merged):,} rows to {OUTPUT_TSV}")


if __name__ == "__main__":
    main()
