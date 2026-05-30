#!/usr/bin/env python3
"""Watch article sources and rebuild only the numbered overlay TeX file."""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import time


def watched_files(root: pathlib.Path, base_pdf: str) -> list[pathlib.Path]:
    article_dir = root / "article"
    files = [article_dir / "main.tex", article_dir / base_pdf]
    stem = pathlib.Path(base_pdf).stem
    files.extend([article_dir / f"{stem}.synctex", article_dir / f"{stem}.synctex.gz"])
    files.extend(sorted((article_dir / "sections").glob("*.tex")))
    return files


def snapshot(paths: list[pathlib.Path]) -> tuple[tuple[str, int, int], ...]:
    state: list[tuple[str, int, int]] = []
    for path in paths:
        try:
            stat = path.stat()
        except FileNotFoundError:
            state.append((path.as_posix(), -1, -1))
            continue
        state.append((path.as_posix(), stat.st_mtime_ns, stat.st_size))
    return tuple(state)


def has_synctex(article_dir: pathlib.Path, base_pdf: str) -> bool:
    stem = pathlib.Path(base_pdf).stem
    return (article_dir / f"{stem}.synctex").exists() or (article_dir / f"{stem}.synctex.gz").exists()


def build_overlay(root: pathlib.Path, base_pdf: str, out_tex: pathlib.Path) -> int:
    article_dir = root / "article"
    if not (article_dir / base_pdf).exists():
        print(f"waiting for article/{base_pdf}", flush=True)
        return 0
    if not has_synctex(article_dir, base_pdf):
        print(f"waiting for SyncTeX file matching article/{base_pdf}", flush=True)
        return 0

    proc = subprocess.run(
        [
            "python3",
            "scripts/build_numbered_overlay.py",
            "--article-dir",
            article_dir.as_posix(),
            "--base-pdf",
            base_pdf,
            "--out-tex",
            out_tex.as_posix(),
        ],
        cwd=root,
        check=False,
    )
    if proc.returncode:
        print(f"numbered overlay build failed with exit code {proc.returncode}", flush=True)
    else:
        print(f"updated {out_tex.relative_to(root).as_posix()}", flush=True)
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[1])
    parser.add_argument("--base-pdf", default="main.pdf")
    parser.add_argument("--out-tex", type=pathlib.Path, default=None)
    parser.add_argument("--interval", type=float, default=1.0)
    parser.add_argument("--debounce", type=float, default=1.5)
    parser.add_argument("--build-on-start", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    article_dir = root / "article"
    out_tex = args.out_tex or article_dir / "main_numbered_overlay.tex"
    out_tex = out_tex.resolve()
    paths = watched_files(root, args.base_pdf)
    last = snapshot(paths)

    print(f"watching {len(paths)} source/PDF files under {root}", flush=True)
    if args.build_on_start:
        build_overlay(root, args.base_pdf, out_tex)
        last = snapshot(paths)

    while True:
        time.sleep(args.interval)
        current = snapshot(paths)
        if current == last:
            continue

        time.sleep(args.debounce)
        target = snapshot(paths)
        print("change detected; rebuilding numbered overlay TeX", flush=True)
        build_overlay(root, args.base_pdf, out_tex)
        last = snapshot(paths)
        if last != target:
            continue


if __name__ == "__main__":
    raise SystemExit(main())
