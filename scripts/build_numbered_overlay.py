#!/usr/bin/env python3
"""Build a LaTeX overlay with physical source-line numbers from SyncTeX."""

from __future__ import annotations

import argparse
import collections
import pathlib
import re
import subprocess


SECTION_INPUT_RE = re.compile(r"\\input\{sections/([^}]+)\}")


def section_files(article_dir: pathlib.Path) -> list[pathlib.Path]:
    main = article_dir / "main.tex"
    files: list[pathlib.Path] = []
    for match in SECTION_INPUT_RE.finditer(main.read_text(encoding="utf-8")):
        files.append(article_dir / "sections" / f"{match.group(1)}.tex")
    return files


def nonblank_lines(path: pathlib.Path) -> list[int]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return [i for i, line in enumerate(lines, start=1) if line.strip()]


def synctex_positions(article_dir: pathlib.Path, pdf_name: str, tex_path: pathlib.Path, line: int) -> list[tuple[int, float, float]]:
    rel = tex_path.relative_to(article_dir).as_posix()
    proc = subprocess.run(
        ["synctex", "view", "-i", f"{line}:1:{rel}", "-o", pdf_name],
        cwd=article_dir,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    positions: list[tuple[int, float, float]] = []
    page: int | None = None
    x: float | None = None
    y: float | None = None
    for raw in proc.stdout.splitlines():
        if raw.startswith("Page:"):
            page = int(raw.split(":", 1)[1])
        elif raw.startswith("x:"):
            x = float(raw.split(":", 1)[1])
        elif raw.startswith("y:"):
            y = float(raw.split(":", 1)[1])
        elif raw == "after:" and page is not None and x is not None and y is not None:
            positions.append((page, x, y))
            page = None
            x = None
            y = None
    return positions


def pick_position(positions: list[tuple[int, float, float]], page_height: float) -> tuple[int, float] | None:
    left_margin_tolerance = 70.0
    filtered = [
        (page, x, y)
        for page, x, y in positions
        if x >= left_margin_tolerance and 80.0 <= y <= page_height - 60.0
    ]
    if not filtered:
        filtered = [
            (page, x, y)
            for page, x, y in positions
            if x >= left_margin_tolerance and 60.0 <= y <= page_height - 40.0
        ]
    if not filtered:
        return None
    page, _x, y = min(filtered, key=lambda item: (item[1], item[0], item[2]))
    return page, y


def compact_numeric_labels(labels: list[str]) -> str:
    nums = sorted(int(label) for label in labels)
    ranges: list[str] = []
    start = nums[0]
    prev = nums[0]
    for num in nums[1:]:
        if num == prev + 1:
            prev = num
            continue
        ranges.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = num
    ranges.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(ranges)


def grouped_marks(marks: list[tuple[float, str]]) -> list[tuple[float, str]]:
    numeric_groups: dict[float, list[str]] = collections.defaultdict(list)
    for y, label in marks:
        numeric_groups[round(y, 1)].append(label)
    grouped = [(y, compact_numeric_labels(labels)) for y, labels in numeric_groups.items()]
    return sorted(grouped, key=lambda item: -item[0])


def pdf_geometry(article_dir: pathlib.Path, base_pdf: str) -> tuple[int, float, float]:
    info = subprocess.run(
        ["pdfinfo", base_pdf],
        cwd=article_dir,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    max_page = 1
    page_width = 612.0
    page_height = 792.0
    for raw in info.stdout.splitlines():
        if raw.startswith("Pages:"):
            max_page = int(raw.split(":", 1)[1].strip())
        elif raw.startswith("Page size:"):
            nums = re.findall(r"[0-9]+(?:\.[0-9]+)?", raw)
            if len(nums) >= 2:
                page_width = float(nums[0])
                page_height = float(nums[1])
    return max_page, page_width, page_height


def build_overlay(article_dir: pathlib.Path, base_pdf: str, out_tex: pathlib.Path) -> None:
    page_marks: dict[int, list[tuple[float, str]]] = collections.defaultdict(list)
    max_page, page_width, page_height = pdf_geometry(article_dir, base_pdf)

    for tex_path in section_files(article_dir):
        for line in nonblank_lines(tex_path):
            pos = pick_position(synctex_positions(article_dir, base_pdf, tex_path, line), page_height)
            if pos is None:
                continue
            page, y_from_top = pos
            page_marks[page].append((page_height - y_from_top, str(line)))
    lines = [
        r"\documentclass{article}",
        rf"\usepackage[paperwidth={page_width:.3f}pt,paperheight={page_height:.3f}pt,margin=0pt]{{geometry}}",
        r"\usepackage{pdfpages}",
        r"\usepackage{xcolor}",
        r"\usepackage{tikz}",
        r"\pagestyle{empty}",
        r"\begin{document}",
    ]

    for page in range(1, max_page + 1):
        marks = grouped_marks(page_marks.get(page, []))
        lines.append(r"\includepdf[pages=%d,pagecommand={},picturecommand={%%" % page)
        lines.append(r"\begin{tikzpicture}[remember picture,overlay]")
        for y, label in marks:
            x = 63.0
            style = r"\tiny\color{red!70!black}"
            lines.append(
                rf"\node[anchor=east] at ([xshift={x:.2f}pt,yshift={y:.2f}pt]current page.south west) {{{style} {label}}};"
            )
        lines.append(r"\end{tikzpicture}%")
        lines.append(r"}]{%s}" % base_pdf)

    lines.append(r"\end{document}")
    out_tex.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--article-dir", type=pathlib.Path, required=True)
    parser.add_argument("--base-pdf", required=True)
    parser.add_argument("--out-tex", type=pathlib.Path, required=True)
    args = parser.parse_args()
    build_overlay(args.article_dir.resolve(), args.base_pdf, args.out_tex.resolve())


if __name__ == "__main__":
    main()
