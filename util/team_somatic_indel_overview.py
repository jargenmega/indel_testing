#!/usr/bin/env python3
"""
Generate overview plots and summary statistics for Team B somatic indels.
Filters the Team B parquet to somatic ALT alleles, normalizes multi-allelic
records into single-indel format, and reports length distributions, chromosome
distributions, and per-base rates for all calls, coding regions, and
not-difficult regions.
"""

from __future__ import annotations

import argparse
import io
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_DATASET = Path(
    "/home/ubuntu/data/teamb_indel-caller/results/Oct15/sample1_indels_only.parquet"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/ubuntu/data/indel_comparison_results/Team_stats"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize Team B somatic indel parquet content."
    )
    parser.add_argument(
        "parquet_path",
        nargs="?",
        type=Path,
        default=DEFAULT_DATASET,
        help=f"Path to Team B indel parquet (default: {DEFAULT_DATASET})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to write summary text and plots (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--exclude-near-germ",
        action="store_true",
        help="Exclude ALT alleles where near_germ == True before generating summaries.",
    )
    return parser.parse_args()


def convert_to_single_indel_format(ref: str, alt: str) -> Tuple[str, str]:
    """Convert REF/ALT pair so exactly one sequence is 1bp (anchor + indel)."""
    anchor = ref[0]
    if len(alt) > len(ref):
        inserted_length = len(alt) - len(ref)
        inserted_bases = alt[1:1 + inserted_length]
        return anchor, anchor + inserted_bases
    if len(alt) < len(ref):
        deleted_length = len(ref) - len(alt)
        deleted_bases = ref[1:1 + deleted_length]
        return anchor + deleted_bases, anchor
    return ref, alt


def format_title(title: str) -> str:
    bar = "=" * len(title)
    return f"{bar}\n{title}\n{bar}"


CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrx": 23, "chry": 24, "chrm": 25, "chrmt": 25})


def chromosome_sort_key(label: str) -> Tuple[int, int, str]:
    """
    Provide a sortable key to ensure chromosomes appear in natural order.
    """
    original = str(label)
    normalized = original.lower()
    if normalized in CHROM_ORDER:
        return (0, CHROM_ORDER[normalized], original)

    stripped = normalized.lstrip("chr")
    if stripped.isdigit():
        return (1, int(stripped), original)

    return (2, float("inf"), original)


def load_dataset(path: Path) -> pd.DataFrame:
    columns = [
        "chrom",
        "pos",
        "allele_seq",
        "allele_type",
        "status",
        "in_coding",
        "in_notdifficult",
        "near_germ",
    ]
    df = pd.read_parquet(path, columns=columns)

    # Map each position to its REF sequence (drop duplicates if any).
    ref_df = (
        df[df["allele_type"] == "REF"]
        .assign(key=lambda d: list(zip(d["chrom"], d["pos"])))
        .drop_duplicates("key")
    )
    ref_map = ref_df.set_index("key")["allele_seq"]

    somatic_alt = df[
        (df["allele_type"] == "ALT") & (df["status"] == "SOMATIC")
    ].copy()
    somatic_alt["key"] = list(zip(somatic_alt["chrom"], somatic_alt["pos"]))
    somatic_alt["ref_seq"] = somatic_alt["key"].map(ref_map)

    missing_ref = somatic_alt["ref_seq"].isna()
    if missing_ref.any():
        missing_positions = (
            somatic_alt.loc[missing_ref, ["chrom", "pos"]]
            .drop_duplicates()
            .to_records(index=False)
        )
        raise ValueError(
            f"Missing REF sequence for positions: {list(missing_positions)}"
        )

    def normalize(row: pd.Series) -> Tuple[str, str]:
        ref_seq = row["ref_seq"]
        alt_seq = row["allele_seq"]
        if len(ref_seq) == 1 or len(alt_seq) == 1:
            return ref_seq, alt_seq
        return convert_to_single_indel_format(ref_seq, alt_seq)

    converted = somatic_alt.apply(normalize, axis=1)
    somatic_alt["converted_ref"] = [pair[0] for pair in converted]
    somatic_alt["converted_alt"] = [pair[1] for pair in converted]

    somatic_alt["ref_len"] = somatic_alt["converted_ref"].str.len()
    somatic_alt["alt_len"] = somatic_alt["converted_alt"].str.len()
    somatic_alt["indel_size"] = (somatic_alt["alt_len"] - somatic_alt["ref_len"]).abs()
    somatic_alt = somatic_alt[somatic_alt["indel_size"] > 0].copy()

    somatic_alt["indel_type"] = np.where(
        somatic_alt["alt_len"] > somatic_alt["ref_len"],
        "insertion",
        "deletion",
    )

    somatic_alt = somatic_alt.drop(columns=["key", "ref_len", "alt_len"])
    return somatic_alt


def indel_length_distribution(df: pd.DataFrame) -> pd.Series:
    counts = (
        df["indel_size"]
        .value_counts(sort=False)
        .sort_index()
        .rename("count")
    )
    return counts


def chromosome_distribution(df: pd.DataFrame) -> pd.Series:
    counts = df["chrom"].value_counts().rename("count")
    return counts.sort_index(key=lambda idx: idx.map(chromosome_sort_key))


def per_base_rates(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "indel_count",
                "unique_positions",
                "span_bp",
                "rate_per_base_total",
            ]
        )

    per_chrom = (
        df.groupby("chrom")
        .agg(
            indel_count=("pos", "size"),
            unique_positions=("pos", "nunique"),
            min_pos=("pos", "min"),
            max_pos=("pos", "max"),
        )
        .reset_index()
    )

    per_chrom = per_chrom.sort_values(
        by="chrom", key=lambda col: col.map(chromosome_sort_key)
    )

    per_chrom["span_bp"] = (
        per_chrom["max_pos"] - per_chrom["min_pos"] + 1
    ).clip(lower=1)

    per_chrom["rate_per_base_total"] = (
        per_chrom["indel_count"] / per_chrom["span_bp"]
    )

    overall = pd.DataFrame(
        {
            "chrom": ["ALL"],
            "indel_count": [per_chrom["indel_count"].sum()],
            "unique_positions": [per_chrom["unique_positions"].sum()],
            "span_bp": [per_chrom["span_bp"].sum()],
        }
    )
    overall["rate_per_base_total"] = (
        overall["indel_count"] / overall["span_bp"]
    )

    rates = pd.concat([per_chrom, overall], ignore_index=True, sort=False)
    return rates[
        [
            "chrom",
            "indel_count",
            "unique_positions",
            "span_bp",
            "rate_per_base_total",
        ]
    ]


class DualWriter:
    """Mirror writes to stdout and an in-memory buffer."""

    def __init__(self) -> None:
        self.buffer = io.StringIO()

    def write(self, text: str = "") -> None:
        print(text)
        self.buffer.write(f"{text}\n")

    def contents(self) -> str:
        return self.buffer.getvalue()


def summarize_subset(name: str, df: pd.DataFrame, writer: DualWriter) -> Tuple[pd.Series, pd.Series, pd.DataFrame]:
    total_calls = len(df)
    unique_positions = df["pos"].nunique()
    writer.write(format_title(f"{name.upper()} SUMMARY"))
    writer.write(f"Total calls: {total_calls:,}")
    writer.write(f"Unique positions: {unique_positions:,}")
    writer.write("")

    length_dist = indel_length_distribution(df)
    writer.write("Length distribution (indel_size -> count):")
    if length_dist.empty:
        writer.write("  <no data>")
    else:
        writer.write(length_dist.to_string())
    writer.write("")

    chrom_dist = chromosome_distribution(df)
    writer.write("Chromosome distribution (chrom -> count):")
    if chrom_dist.empty:
        writer.write("  <no data>")
    else:
        writer.write(chrom_dist.to_string())
    writer.write("")

    rates = per_base_rates(df)
    writer.write("Per-base rates:")
    if rates.empty:
        writer.write("  <no data>")
    else:
        writer.write(rates.to_string(index=False))
    writer.write("")
    writer.write("-" * 80)
    writer.write("")

    return length_dist, chrom_dist, rates


def apply_log_if_appropriate(ax, values: Iterable[float], axis: str = "y") -> None:
    values = np.asarray(list(values), dtype=float)
    positive = values[values > 0]
    if positive.size == 0:
        return
    ratio = positive.max() / positive.min() if positive.min() > 0 else np.inf
    if ratio >= 100:
        getattr(ax, f"set_{axis}scale")("log")


def safe_name(name: str) -> str:
    return name.lower().replace(" ", "_")


def plot_length_distribution(name: str, counts: pd.Series, output_dir: Path) -> Path:
    if counts.empty:
        return None

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(counts.index, counts.values, drawstyle="steps-mid")
    ax.set_title(f"Indel Length Distribution ({name})")
    ax.set_xlabel("Indel size (bp)")
    ax.set_ylabel("Count")
    apply_log_if_appropriate(ax, counts.values, axis="y")
    fig.tight_layout()
    filename = output_dir / f"{safe_name(name)}_length_distribution.png"
    fig.savefig(filename)
    plt.close(fig)
    return filename


def plot_chromosome_distribution(name: str, counts: pd.Series, output_dir: Path) -> Path:
    if counts.empty:
        return None

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(counts.index, counts.values)
    ax.set_title(f"Chromosome Distribution ({name})")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Count")
    ax.tick_params(axis="x", rotation=45)
    for tick in ax.get_xticklabels():
        tick.set_horizontalalignment("right")
    apply_log_if_appropriate(ax, counts.values, axis="y")
    fig.tight_layout()
    filename = output_dir / f"{safe_name(name)}_chromosome_distribution.png"
    fig.savefig(filename)
    plt.close(fig)
    return filename


def plot_per_base_rates(name: str, rates: pd.DataFrame, output_dir: Path) -> Path:
    data = rates[rates["chrom"] != "ALL"].copy()
    if data.empty:
        return None

    data = data.sort_values(
        by="chrom", key=lambda col: col.map(chromosome_sort_key)
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(
        data["chrom"],
        data["rate_per_base_total"],
        marker="o",
    )
    ax.set_title(f"Per-base Indel Rate ({name})")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Rate per base")
    ax.tick_params(axis="x", rotation=45)
    for tick in ax.get_xticklabels():
        tick.set_horizontalalignment("right")
    apply_log_if_appropriate(
        ax,
        data["rate_per_base_total"].values,
        axis="y",
    )
    fig.tight_layout()
    filename = output_dir / f"{safe_name(name)}_per_base_rate.png"
    fig.savefig(filename)
    plt.close(fig)
    return filename


def main() -> None:
    args = parse_args()
    parquet_path = args.parquet_path
    output_dir = args.output_dir

    if not parquet_path.exists():
        raise FileNotFoundError(f"Parquet file not found: {parquet_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    writer = DualWriter()

    writer.write(f"Loading {parquet_path} ...")
    df = load_dataset(parquet_path)
    if args.exclude_near_germ:
        df = df[df["near_germ"] == False].copy()  # noqa: E712
        writer.write("Applied filter: near_germ == False")
    writer.write(f"Loaded {len(df):,} somatic ALT rows")
    writer.write("")

    subsets: Dict[str, pd.DataFrame] = {
        "all positions": df,
        "coding": df[df["in_coding"]],
        "not difficult": df[df["in_notdifficult"]],
    }

    generated_charts: List[Path] = []

    for name, subset in subsets.items():
        length_dist, chrom_dist, rates = summarize_subset(name, subset, writer)
        chart = plot_length_distribution(name, length_dist, output_dir)
        if chart:
            generated_charts.append(chart)
        chart = plot_chromosome_distribution(name, chrom_dist, output_dir)
        if chart:
            generated_charts.append(chart)
        chart = plot_per_base_rates(name, rates, output_dir)
        if chart:
            generated_charts.append(chart)

    summary_path = output_dir / f"{parquet_path.stem}_somatic_summary.txt"
    summary_path.write_text(writer.contents())

    print(f"Summary written to {summary_path}")
    if generated_charts:
        print("Generated charts:")
        for chart_path in generated_charts:
            print(f"  {chart_path}")
    else:
        print("No charts generated.")


if __name__ == "__main__":
    main()
