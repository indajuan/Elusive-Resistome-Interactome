#!/usr/bin/env python3
"""
Port of the R core-resistome analysis:

  cut_size_core / filter_samples_core / core_resistome / the seeds loop

For `number_of_subsamples` iterations, this repeatedly draws
`sub_sampling_size` metagenome samples per habitat (without replacement),
and for each ARG-detection tool asks: for each vsearch cluster centroid,
in what proportion (p) of the drawn samples for a habitat does it appear?
Centroids with p >= cut (for cut in 0.2..0.9) are tallied. After all
iterations, "cnt" is the number of iterations (out of number_of_subsamples)
in which a (centroid, new_level_centroid, tool, habitat) passed each cut --
i.e. how consistently it qualifies as "core" resistome at that threshold.

Differences from the R code (by design, see conversation):
  - habitats are restricted to the 13 used for the unigenes habitat filter
    (drops amplicon / built-environment / isolate from metadata).
  - habitats with fewer than sub_sampling_size samples (e.g. wastewater)
    draw all available samples instead of erroring
    (n = min(sub_sampling_size, samples available in that habitat)).
  - the per-iteration seed is seed, seed+1, seed+2, ... (sequential from
    the same seed passed to rarefaction), rather than an externally
    supplied seed list.
  - the R code's per-iteration bind_rows + re-summarise(cnt=sum(cnt)) is
    algebraically equivalent to accumulating every iteration's hits and
    summing once at the end; this does the latter for speed.

Usage:
  ./core_resistome.py [data_dir] [output_tsv] [sub_sampling_size] [number_of_subsamples] [depth] [seed]

  data_dir              defaults to data_zenodo_github
  output_tsv            defaults to core_resistome.tsv
  sub_sampling_size     defaults to 100
  number_of_subsamples  defaults to 500
  depth                 defaults to 5e6   (rarefaction depth)
  seed                  defaults to 2000  (rarefaction seed; also seeds
                                            the per-iteration subsampling
                                            sequentially: seed, seed+1, ...)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from rarefy_abundances import rarefy

CUTS = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

HABITATS = [
    "human gut", "human oral", "human skin", "human nose", "human vagina",
    "dog gut", "cat gut", "mouse gut", "pig gut",
    "wastewater", "marine", "freshwater", "soil",
]


def load_clusters(data_dir: Path) -> dict:
    """clusters.uc (vsearch UC format) -> {query_id: centroid_id}."""
    cols = [f"V{i}" for i in range(1, 11)]
    df = pd.read_csv(data_dir / "clusters.uc", sep="\t", header=None, names=cols, dtype=str)
    df = df[df["V1"] != "C"]
    df["centroid"] = df["V10"].where(df["V10"] != "*", df["V9"])
    return dict(zip(df["V9"], df["centroid"]))


def build_lst(unigenes: pd.DataFrame) -> pd.DataFrame:
    """unigenes.tsv, deduplicated to one new_level_centroid per (tool, query) --
    mirrors R's match(), which keeps the first hit."""
    lst = unigenes[["query", "tool", "new_level_centroid"]].drop_duplicates(
        subset=["tool", "query"], keep="first"
    )
    return lst


def prepare_args_abundances(data_dir: Path, metadata: pd.DataFrame, centroid_map: dict,
                             depth: float, seed: int) -> pd.DataFrame:
    print("Reading args_abundances.tsv.gz ...")
    args_abundances = pd.read_csv(data_dir / "args_abundances.tsv.gz", sep="\t")
    if "X" not in args_abundances.columns and "Unnamed: 0" in args_abundances.columns:
        # the file's first column (gene ID) has a blank header, so pandas
        # names it "Unnamed: 0" instead of "X" -- recover the real name
        args_abundances = args_abundances.rename(columns={"Unnamed: 0": "X"})

    print(f"Rarefying to depth={depth:,.0f}, seed={seed} ...")
    args_abundances = rarefy(args_abundances, metadata, depth=depth, seed=seed)
    args_abundances = args_abundances[args_abundances["rarified_count"] > 0]
    args_abundances = args_abundances.rename(columns={"X": "query"})

    args_abundances["centroid"] = args_abundances["query"].map(centroid_map)
    args_abundances = args_abundances.dropna(subset=["centroid"])

    habitat_map = metadata.set_index("sample_id")["habitat"].to_dict()
    args_abundances["habitat"] = args_abundances["sample"].map(habitat_map)
    args_abundances = args_abundances[args_abundances["habitat"].isin(HABITATS)]

    return args_abundances[["query", "sample", "centroid", "habitat"]]


def draw_samples(rng: np.random.Generator, samples_to_collect: pd.DataFrame, sub_sampling_size: int):
    """Per habitat, draw min(sub_sampling_size, available) samples without
    replacement. Returns (list of drawn sample ids, {habitat: N})."""
    drawn = []
    counts = {}
    for habitat, group in samples_to_collect.groupby("habitat"):
        available = group["sample"].to_numpy()
        n = min(sub_sampling_size, len(available))
        chosen = rng.choice(available, size=n, replace=False)
        drawn.append(chosen)
        counts[habitat] = n
    return np.concatenate(drawn), counts


def core_resistome_iteration(args_abundances: pd.DataFrame, lst: pd.DataFrame,
                              rng: np.random.Generator, samples_to_collect: pd.DataFrame,
                              sub_sampling_size: int, cuts) -> pd.DataFrame:
    drawn_samples, samples_to_collect_number = draw_samples(rng, samples_to_collect, sub_sampling_size)

    args_abundances_core = args_abundances[args_abundances["sample"].isin(set(drawn_samples))]

    # one merge instead of looping per tool (equivalent to rbind of filter_samples_core over lst)
    hits = args_abundances_core.merge(lst, on="query", how="inner")
    assert hits["new_level_centroid"].notna().all(), "NA values detected in new_level_centroid"

    grp = hits.groupby(["tool", "habitat", "centroid"], as_index=False).agg(
        new_level_centroid=("new_level_centroid", "first"),
        n=("sample", "nunique"),
    )
    grp["N"] = grp["habitat"].map(samples_to_collect_number)
    grp["p"] = grp["n"] / grp["N"]
    grp = grp[grp["p"] >= 0.1]

    cuts_arr = np.asarray(cuts)
    satisfied = grp["p"].to_numpy()[:, None] >= cuts_arr[None, :]
    rows_idx, cuts_idx = np.nonzero(satisfied)

    return pd.DataFrame({
        "centroid": grp["centroid"].to_numpy()[rows_idx],
        "new_level_centroid": grp["new_level_centroid"].to_numpy()[rows_idx],
        "tool": grp["tool"].to_numpy()[rows_idx],
        "habitat": grp["habitat"].to_numpy()[rows_idx],
        "cut": cuts_arr[cuts_idx],
        "cnt": 1,
    })


def run_core_resistome(data_dir: Path, sub_sampling_size: int, number_of_subsamples: int,
                        depth: float = 5e6, seed: int = 2000, cuts=CUTS) -> pd.DataFrame:
    metadata = pd.read_csv(data_dir / "metagenomes_metadata.csv")
    metadata = metadata[metadata["habitat"].isin(HABITATS)]

    print("Loading clusters.uc ...")
    centroid_map = load_clusters(data_dir)

    print("Loading unigenes.tsv ...")
    unigenes = pd.read_csv(data_dir.parent / "unigenes.tsv", sep="\t")
    lst = build_lst(unigenes)

    args_abundances = prepare_args_abundances(data_dir, metadata, centroid_map, depth, seed)
    print(f"  {len(args_abundances):,} rarefied, centroid-mapped, habitat-restricted abundance rows")

    samples_to_collect = metadata[["sample_id", "habitat"]].drop_duplicates().rename(columns={"sample_id": "sample"})

    iterations = []
    for j in range(number_of_subsamples):
        sed = seed + j
        print(f"Subsample iteration {j + 1}/{number_of_subsamples} (seed={sed})")
        rng = np.random.default_rng(sed)
        iterations.append(
            core_resistome_iteration(args_abundances, lst, rng, samples_to_collect, sub_sampling_size, cuts)
        )

    combined = pd.concat(iterations, ignore_index=True)
    result = combined.groupby(
        ["centroid", "new_level_centroid", "tool", "habitat", "cut"], as_index=False
    )["cnt"].sum()
    return result


def main():
    data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data_zenodo_github")
    output_tsv = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("core_resistome.tsv")
    sub_sampling_size = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    number_of_subsamples = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    depth = float(sys.argv[5]) if len(sys.argv) > 5 else 5e6
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 2000

    result = run_core_resistome(data_dir, sub_sampling_size, number_of_subsamples, depth, seed)

    result.to_csv(output_tsv, sep="\t", index=False)
    print(f"Wrote {len(result):,} rows to {output_tsv}")


if __name__ == "__main__":
    main()
