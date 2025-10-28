#!/usr/bin/env python3
"""
Count Team B indel calls that land immediately downstream of difficult regions.

For each ALT allele (default: somatic) in the provided Team B parquet, the script
checks whether the indel position lies +1 to +N bases after the end of any interval
in the GRCh38 difficult-region union BED. It reports counts for each offset bucket,
prints summary totals, and optionally plots the distribution.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_PARQUET = Path(
    "/home/ubuntu/data/teamb_indel-caller/results/Oct15/sample1_indels_only.parquet"
)
DEFAULT_BED = Path(
    "/home/ubuntu/data/reference/beds/GRCh38_difficult_union.bed.gz"
)
DEFAULT_WINDOW = 50


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count Team B indels immediately downstream of difficult-region ends."
    )
    parser.add_argument(
        "--team-parquet",
        type=Path,
        default=DEFAULT_PARQUET,
        help=f"Team B parquet with indel alleles (default: {DEFAULT_PARQUET})",
    )
    parser.add_argument(
        "--bed-path",
        type=Path,
        default=DEFAULT_BED,
        help=f"Difficult-region BED (default: {DEFAULT_BED})",
    )
    parser.add_argument(
        "--status-filter",
        type=str,
        default="SOMATIC",
        choices=["SOMATIC", "GERMLINE", "ALL"],
        help="Status to include for ALT alleles (default: SOMATIC).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional path to write the summary table (TSV).",
    )
    parser.add_argument(
        "--plot-output",
        type=Path,
        help="Optional path to write a plot (PNG) of distance vs count.",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        help=f"Maximum distance (bp) downstream to analyze (default: {DEFAULT_WINDOW}).",
    )
    parser.add_argument(
        "--require-not-near-germ",
        action="store_true",
        help="Only include ALT alleles where near_germ == False.",
    )
    parser.add_argument(
        "--require-in-notdifficult",
        action="store_true",
        help="Only include ALT alleles where in_notdifficult == True.",
    )
    return parser.parse_args()


def load_bed_intervals(bed_path: Path) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Load BED intervals keyed by chromosome.
    Returns dict: chrom -> {'start0': np.ndarray, 'end0': np.ndarray}
    where coordinates remain 0-based half-open.
    """
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    bed_df = pd.read_csv(
        bed_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start0", "end0"],
    )

    intervals: Dict[str, Dict[str, np.ndarray]] = {}
    for chrom, group in bed_df.groupby("chrom"):
        sorted_group = group.sort_values("start0").reset_index(drop=True)
        intervals[chrom] = {
            "start0": sorted_group["start0"].to_numpy(dtype=np.int64),
            "end0": sorted_group["end0"].to_numpy(dtype=np.int64),
        }
    return intervals


def load_team_alt_alleles(
    parquet_path: Path,
    status_filter: str,
    require_not_near_germ: bool,
    require_in_notdifficult: bool,
) -> pd.DataFrame:
    """
    Load ALT alleles from Team B parquet with optional status filtering.
    """
    if not parquet_path.exists():
        raise FileNotFoundError(f"Team parquet not found: {parquet_path}")

    columns = [
        "chrom",
        "pos",
        "allele_type",
        "status",
        "near_germ",
        "in_notdifficult",
    ]
    df = pd.read_parquet(parquet_path, columns=columns)

    alt_df = df[df["allele_type"] == "ALT"].copy()
    if status_filter != "ALL":
        alt_df = alt_df[alt_df["status"] == status_filter].copy()

    if require_not_near_germ and "near_germ" in alt_df.columns:
        alt_df = alt_df[alt_df["near_germ"] == False].copy()  # noqa: E712

    if require_in_notdifficult and "in_notdifficult" in alt_df.columns:
        alt_df = alt_df[alt_df["in_notdifficult"] == True].copy()  # noqa: E712

    alt_df["pos"] = alt_df["pos"].astype(np.int64)
    return alt_df


def count_offsets(
    alt_df: pd.DataFrame,
    intervals: Dict[str, Dict[str, np.ndarray]],
    offsets: Iterable[int],
) -> Dict[int, int]:
    """
    Count how many alt alleles fall within each offset (1..N) downstream of
    a difficult-region end.
    """
    counts = {offset: 0 for offset in offsets}

    for chrom, chrom_df in alt_df.groupby("chrom"):
        if chrom not in intervals:
            continue

        starts0 = intervals[chrom]["start0"]
        ends0 = intervals[chrom]["end0"]
        if starts0.size == 0:
            continue

        positions = chrom_df["pos"].to_numpy(dtype=np.int64)
        pos0 = positions - 1  # convert to 0-based coordinates

        # Find the first interval start that is strictly greater than pos0
        idx_right = np.searchsorted(starts0, pos0, side="right")

        for pos_1b, pos0_val, idx in zip(positions, pos0, idx_right):
            # Determine if the position is still inside an interval
            in_interval = idx > 0 and pos0_val < ends0[idx - 1]
            if in_interval:
                continue  # inside a difficult region, skip

            # We care about the previous interval (the one that just ended before this position)
            if idx == 0:
                continue  # there is no interval ending before this position

            end0 = ends0[idx - 1]  # end coordinate is exclusive (0-based)
            end_1b = end0  # convert to 1-based
            if pos_1b <= end_1b:
                continue  # inside or before the end boundary

            diff = pos_1b - end_1b
            if diff in counts:
                counts[diff] += 1

    counts["total_alt_evaluated"] = len(alt_df)
    counts["total_near_boundary"] = sum(counts[offset] for offset in offsets)
    return counts


def main() -> None:
    args = parse_args()

    intervals = load_bed_intervals(args.bed_path)
    alt_df = load_team_alt_alleles(
        args.team_parquet,
        args.status_filter,
        args.require_not_near_germ,
        args.require_in_notdifficult,
    )

    if alt_df.empty:
        print("No ALT alleles match the requested filters.")
        return

    offsets = range(1, args.window + 1)
    counts = count_offsets(alt_df, intervals, offsets)

    summary_rows = []
    for offset in offsets:
        summary_rows.append(
            {
                "distance_bp_from_end": offset,
                "count": counts[offset],
            }
        )

    summary_df = pd.DataFrame(summary_rows)

    print("\nCounts of ALT alleles downstream of difficult-region ends")
    print(f"Team parquet: {args.team_parquet}")
    print(f"Status filter: {args.status_filter}")
    print(summary_df.to_string(index=False))

    total_alt = counts["total_alt_evaluated"]
    total_near = counts["total_near_boundary"]
    pct_near = (total_near / total_alt * 100) if total_alt else 0.0

    print(f"\nTotal ALT alleles evaluated: {total_alt}")
    print(f"Total within 1-{args.window}bp of a difficult-region end: {total_near}")
    print(f"Percentage within window: {pct_near:.2f}%")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w") as handle:
            summary_df.to_csv(handle, sep="\t", index=False)
            handle.write("\n")
            handle.write("# Filters applied:\n")
            handle.write(f"#   status_filter = {args.status_filter}\n")
            handle.write(f"#   near_germ == False? {args.require_not_near_germ}\n")
            handle.write(f"#   in_notdifficult == True? {args.require_in_notdifficult}\n")
            handle.write(f"# Window size (bp): {args.window}\n")
            handle.write(
                f"# Total ALT alleles evaluated: {total_alt}\n"
            )
            handle.write(
                f"# Total within {args.window}bp: {total_near} "
                f"({pct_near:.2f}% of evaluated)\n"
            )
        print(f"\nSummary written to {args.output}")

    plot_path = args.plot_output
    if plot_path is None and args.output is not None:
        plot_path = args.output.with_suffix(".png")

    if plot_path:
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        distances = summary_df["distance_bp_from_end"]
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.bar(distances, summary_df["count"], color="tab:blue")
        ax.set_xlabel("Distance from difficult-region end (bp)")
        ax.set_ylabel("Indel count")
        ax.set_title("Somatic ALT indels downstream of difficult-region ends")
        major_xticks = distances[distances % 5 == 0]
        ax.set_xticks(major_xticks)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(plot_path)
        plt.close(fig)
        print(f"Plot written to {plot_path}")


if __name__ == "__main__":
    main()
