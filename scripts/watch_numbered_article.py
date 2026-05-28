#!/usr/bin/env python3
"""Watch article TeX sources and rebuild the numbered review PDF."""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import time


def watched_files(root: pathlib.Path) -> list[pathlib.Path]:
    article_dir = root / "article"
    files = [article_dir / "main.tex"]
    files.extend(sorted((article_dir / "sections").glob("*.tex")))
    return [path for path in files if path.exists()]


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[1])
    parser.add_argument("--interval", type=float, default=1.0)
    parser.add_argument("--debounce", type=float, default=1.5)
    parser.add_argument("--compile-on-start", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    compile_script = root / "scripts" / "compile_numbered_article.sh"
    paths = watched_files(root)
    last = snapshot(paths)

    print(f"watching {len(paths)} TeX files under {root}", flush=True)
    if args.compile_on_start:
        subprocess.run([compile_script.as_posix()], cwd=root, check=False)
        last = snapshot(paths)

    while True:
        time.sleep(args.interval)
        current = snapshot(paths)
        if current == last:
            continue

        time.sleep(args.debounce)
        target = snapshot(paths)
        print("change detected; rebuilding article/main_numbered.pdf", flush=True)
        proc = subprocess.run([compile_script.as_posix()], cwd=root, check=False)
        if proc.returncode:
            print(f"numbered article build failed with exit code {proc.returncode}", flush=True)
        after = snapshot(paths)
        last = target
        if after != target:
            continue


if __name__ == "__main__":
    raise SystemExit(main())
