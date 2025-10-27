#!/usr/bin/env python3
"""
Generate summary statistics for a DeepVariant indel parquet file.
Focuses on indel length distribution, chromosome distribution, and per-base rates
for all calls as well as calls in coding and not-difficult regions.
"""

from __future__ import annotations

import argparse
import io
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

DEFAULT_DATASET = Path(
    "/home/ubuntu/data/teamb_indel-caller/deep_variant/GTEx-sample1.indels.parquet"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/ubuntu/data/indel_comparison_results/DV_stats"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize DeepVariant indel parquet content."
    )
    parser.add_argument(
        "parquet_path",
        nargs="?",
        type=Path,
        default=DEFAULT_DATASET,
        help=f"Path to DeepVariant indel parquet (default: {DEFAULT_DATASET})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to write summary text and plots (default: {DEFAULT_OUTPUT_DIR})",
    )
    return parser.parse_args()


def load_dataset(path: Path) -> pd.DataFrame:
    columns = [
        "chrom",
        "pos",
        "indel_size",
        "indel_type",
        "in_coding",
        "in_notdifficult",
    ]
    df = pd.read_parquet(path, columns=columns)
    # Defensive copy to avoid pandas SettingWithCopy warnings downstream.
    return df.copy()


def format_title(title: str) -> str:
    bar = "=" * len(title)
    return f"{bar}\n{title}\n{bar}"


def indel_length_distribution(df: pd.DataFrame) -> pd.Series:
    counts = (
        df["indel_size"]
        .value_counts(sort=False)
        .sort_index()
        .rename("count")
    )
    return counts


def chromosome_distribution(df: pd.DataFrame) -> pd.Series:
    return (
        df["chrom"]
        .value_counts()
        .sort_index()
        .rename("count")
    )


def per_base_rates(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "indel_count",
                "unique_positions",
                "span_bp",
                "rate_per_base_total",
                "rate_per_base_unique",
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

    per_chrom["span_bp"] = (
        per_chrom["max_pos"] - per_chrom["min_pos"] + 1
    ).clip(lower=1)

    per_chrom["rate_per_base_total"] = (
        per_chrom["indel_count"] / per_chrom["span_bp"]
    )
    per_chrom["rate_per_base_unique"] = (
        per_chrom["unique_positions"] / per_chrom["span_bp"]
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
    overall["rate_per_base_unique"] = (
        overall["unique_positions"] / overall["span_bp"]
    )

    rates = pd.concat([per_chrom, overall], ignore_index=True, sort=False)
    return rates[
        [
            "chrom",
            "indel_count",
            "unique_positions",
            "span_bp",
            "rate_per_base_total",
            "rate_per_base_unique",
        ]
    ]


def print_section(title: str, items: Iterable[str]) -> None:
    print(format_title(title))
    for item in items:
        print(item)
    print()

class DualWriter:
    """Mirror writes to stdout and in-memory buffer for later file output."""

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
    indel_types = sorted(df["indel_type"].unique()) if total_calls else []
    writer.write(f"Indel types: {indel_types}")
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

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(
        data["chrom"],
        data["rate_per_base_total"],
        marker="o",
        label="Count / span",
    )
    ax.plot(
        data["chrom"],
        data["rate_per_base_unique"],
        marker="s",
        label="Unique positions / span",
    )
    ax.set_title(f"Per-base Indel Rate ({name})")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Rate per base")
    ax.tick_params(axis="x", rotation=45)
    for tick in ax.get_xticklabels():
        tick.set_horizontalalignment("right")
    apply_log_if_appropriate(
        ax,
        np.concatenate(
            [data["rate_per_base_total"].values, data["rate_per_base_unique"].values]
        ),
        axis="y",
    )
    ax.legend()
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
    writer.write(f"Loaded {len(df):,} rows")
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

    summary_path = output_dir / f"{parquet_path.stem}_summary.txt"
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
