#!/usr/bin/env python3
r"""Generate and validate larger DDGP rank-count instances.

This runner creates moderately pruned DDGP instances up to 30 vertices and then
checks the global labelled-violation rank-count prediction against the full
Branch-and-Prune count from scripts/ddgp_bp.py.

The pruning counts are intentionally high enough to keep exact BP validation
practical. The rank-count computation itself is cheap; exact BP is the limiting
step.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SuiteCase:
    k: int
    n: int
    pruning_count: int
    seed: int

    @property
    def filename(self) -> str:
        return f"ddgp-k{self.k}-n{self.n}-p{self.pruning_count}.jsonl"


DEFAULT_CASES = (
    SuiteCase(k=2, n=20, pruning_count=8, seed=3020),
    SuiteCase(k=2, n=25, pruning_count=12, seed=3025),
    SuiteCase(k=2, n=30, pruning_count=16, seed=3030),
    SuiteCase(k=3, n=20, pruning_count=8, seed=3320),
    SuiteCase(k=3, n=25, pruning_count=12, seed=3325),
    SuiteCase(k=3, n=30, pruning_count=16, seed=3330),
)


def run(command: list[str]) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("research/large-rank-count-validation"),
        help="Directory for generated JSONL datasets.",
    )
    parser.add_argument(
        "--samples-per-case",
        type=int,
        default=5,
        help="Number of generated instances per (K, n, pruning-count) case.",
    )
    parser.add_argument(
        "--count",
        type=int,
        default=None,
        help="Maximum instances to validate per generated file. Defaults to all samples.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-5,
        help="BP pruning tolerance passed to verify_compatibility.py.",
    )
    parser.add_argument(
        "--skip-generate",
        action="store_true",
        help="Reuse existing JSONL files in --out-dir.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    validate_count = args.samples_per_case if args.count is None else args.count

    for case in DEFAULT_CASES:
        path = args.out_dir / case.filename
        if not args.skip_generate:
            run(
                [
                    sys.executable,
                    "scripts/ddgp_generate_multi_pruning.py",
                    "--out",
                    str(path),
                    "--n",
                    str(case.n),
                    "--k",
                    str(case.k),
                    "--mode",
                    "ddgp",
                    "--samples-per-count",
                    str(args.samples_per_case),
                    "--pruning-counts",
                    str(case.pruning_count),
                    "--seed",
                    str(case.seed),
                    "--global-min-area2",
                    "1e-8",
                ]
            )

        run(
            [
                sys.executable,
                "scripts/verify_compatibility.py",
                "--instances",
                str(path),
                "--count",
                str(validate_count),
                "--tol",
                str(args.tol),
            ]
        )


if __name__ == "__main__":
    main()
