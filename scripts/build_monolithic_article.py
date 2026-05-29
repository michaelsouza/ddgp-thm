#!/usr/bin/env python3
"""Build a single-file LaTeX version of article/main.tex."""

from __future__ import annotations

import argparse
import pathlib
import re


INPUT_RE = re.compile(
    r"^(?P<indent>[ \t]*)\\(?P<command>input|include)\s*\{(?P<path>[^}]+)\}\s*(?P<comment>%.*)?$"
)
BIBLIOGRAPHY_RE = re.compile(
    r"^(?P<indent>[ \t]*)\\bibliography\s*\{(?P<path>[^}]+)\}\s*(?P<comment>%.*)?$"
)


class MonolithicBuildError(RuntimeError):
    """Raised when an included LaTeX file cannot be expanded."""


def resolve_tex_path(parent: pathlib.Path, raw_path: str) -> pathlib.Path:
    path = pathlib.Path(raw_path)
    if not path.suffix:
        path = path.with_suffix(".tex")
    if not path.is_absolute():
        path = parent / path
    return path.resolve()


def display_path(path: pathlib.Path, base_dir: pathlib.Path) -> str:
    try:
        return path.relative_to(base_dir).as_posix()
    except ValueError:
        return path.as_posix()


def comment_block(title: str, path: pathlib.Path, base_dir: pathlib.Path) -> list[str]:
    return [
        "",
        f"% === BEGIN {title}: {display_path(path, base_dir)} ===",
    ]


def end_comment_block(title: str, path: pathlib.Path, base_dir: pathlib.Path) -> list[str]:
    return [
        f"% === END {title}: {display_path(path, base_dir)} ===",
        "",
    ]


def expand_tex_file(
    path: pathlib.Path,
    *,
    base_dir: pathlib.Path,
    annotate: bool = False,
    stack: tuple[pathlib.Path, ...] = (),
) -> list[str]:
    path = path.resolve()
    if path in stack:
        cycle = " -> ".join(p.as_posix() for p in (*stack, path))
        raise MonolithicBuildError(f"recursive LaTeX include detected: {cycle}")
    if not path.is_file():
        raise MonolithicBuildError(f"included file not found: {path}")

    lines: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        match = INPUT_RE.match(raw_line)
        if match is None:
            lines.append(raw_line)
            continue

        included = resolve_tex_path(path.parent, match.group("path"))
        if annotate:
            lines.extend(comment_block("INPUT", included, base_dir))
        lines.extend(expand_tex_file(included, base_dir=base_dir, annotate=annotate, stack=(*stack, path)))
        if annotate:
            lines.extend(end_comment_block("INPUT", included, base_dir))

    return lines


def inline_bibliography(
    lines: list[str],
    main_path: pathlib.Path,
    bbl_path: pathlib.Path | None,
    *,
    annotate: bool = False,
) -> list[str]:
    if bbl_path is None:
        bbl_path = main_path.with_suffix(".bbl")
    bbl_path = bbl_path.resolve()
    if not bbl_path.is_file():
        raise MonolithicBuildError(f"bibliography file not found: {bbl_path}")

    bbl_lines = bbl_path.read_text(encoding="utf-8").splitlines()
    out: list[str] = []
    for line in lines:
        if BIBLIOGRAPHY_RE.match(line):
            if annotate:
                out.extend(comment_block("BIBLIOGRAPHY", bbl_path, main_path.parent))
            out.extend(bbl_lines)
            if annotate:
                out.extend(end_comment_block("BIBLIOGRAPHY", bbl_path, main_path.parent))
        else:
            out.append(line)
    return out


def build_monolithic(
    main_path: pathlib.Path,
    out_path: pathlib.Path,
    *,
    inline_bbl: bool = False,
    bbl_path: pathlib.Path | None = None,
    annotate: bool = False,
) -> None:
    main_path = main_path.resolve()
    out_path = out_path.resolve()
    lines = expand_tex_file(main_path, base_dir=main_path.parent, annotate=annotate)
    if inline_bbl:
        lines = inline_bibliography(lines, main_path, bbl_path, annotate=annotate)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Expand article/main.tex into a single LaTeX source file."
    )
    parser.add_argument(
        "--main",
        type=pathlib.Path,
        default=pathlib.Path("article/main.tex"),
        help="Main LaTeX file to expand. Default: article/main.tex",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("article/main_monolithic.tex"),
        help="Output path for the monolithic LaTeX file. Default: article/main_monolithic.tex",
    )
    parser.add_argument(
        "--inline-bbl",
        action="store_true",
        help="Replace \\bibliography{...} with the generated .bbl contents.",
    )
    parser.add_argument(
        "--bbl",
        type=pathlib.Path,
        default=None,
        help="Bibliography .bbl path used with --inline-bbl. Default: sibling .bbl next to --main.",
    )
    parser.add_argument(
        "--annotate",
        action="store_true",
        help="Add comments marking expanded inputs and inlined bibliography using relative paths.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        build_monolithic(
            args.main,
            args.output,
            inline_bbl=args.inline_bbl,
            bbl_path=args.bbl,
            annotate=args.annotate,
        )
    except MonolithicBuildError as exc:
        raise SystemExit(f"error: {exc}") from exc
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
