#!/usr/bin/env python3
"""Compact Scopus CSV exports into deduplicated screening batches."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEARCH_DIR = ROOT / "references" / "searches"
SCREEN_DIR = ROOT / "references" / "screening"
BATCH_DIR = SCREEN_DIR / "batches"
RESULT_DIR = SCREEN_DIR / "agent-results"
BATCH_SIZE = 20
ABSTRACT_LIMIT = 1100

TAG_RULES = [
    ("DDGP/DMDGP definitions", [
        "ddgp",
        "dmdgp",
        "discretizable distance geometry",
        "discretizable molecular distance geometry",
    ]),
    ("Branch-and-Prune", [
        "branch-and-prune",
        "branch and prune",
        "branch-and-bound",
    ]),
    ("partial reflection symmetries", [
        "partial reflection",
        "reflection symmetr",
        "symmetry",
        "symmetries",
        "symmetry-based",
    ]),
    ("solution counting", [
        "solution count",
        "number of solutions",
        "power of two",
        "counting method",
        "combinatorial counting",
    ]),
    ("distance geometry applications", [
        "distance geometry",
        "molecular",
        "protein",
        "sensor",
        "localization",
    ]),
    ("algebraic/genericity methods", [
        "generic",
        "genericity",
        "algebraic",
        "cayley-menger",
        "clifford",
        "geometric algebra",
        "rank",
    ]),
]

KEYWORD_WEIGHTS = {
    "discretizable distance geometry problem": 8,
    "discretizable molecular distance geometry problem": 7,
    "ddgp": 7,
    "dmdgp": 6,
    "branch-and-prune": 5,
    "branch and prune": 5,
    "partial reflection": 5,
    "symmetry": 3,
    "symmetries": 3,
    "solution count": 5,
    "number of solutions": 5,
    "power of two": 4,
    "combinatorial counting": 4,
    "lateration": 2,
    "sphere intersection": 2,
    "distance geometry": 2,
}


def clean_space(value: str | None) -> str:
    return re.sub(r"\s+", " ", value or "").strip()


def normalize_title(value: str | None) -> str:
    value = (value or "").casefold()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return clean_space(value)


def query_id(path: Path) -> str:
    name = path.name
    name = name.removeprefix("scopus-")
    name = name.removesuffix(".csv")
    name = name.removesuffix(".csv")
    name = re.sub(r"-2026-.*$", "", name)
    return name


def truncate(value: str, limit: int = ABSTRACT_LIMIT) -> str:
    value = clean_space(value)
    if len(value) <= limit:
        return value
    cut = value[:limit].rsplit(" ", 1)[0]
    return cut + " ..."


def tags_for(record: dict[str, str]) -> list[str]:
    text = " ".join([
        record.get("Title", ""),
        record.get("Abstract", ""),
        record.get("Author Keywords", ""),
        record.get("Index Keywords", ""),
        record.get("Source title", ""),
    ]).casefold()
    tags = []
    for tag, needles in TAG_RULES:
        if any(needle in text for needle in needles):
            tags.append(tag)
    return tags


def score_record(record: dict[str, str], query_hits: set[str]) -> int:
    text = " ".join([
        record.get("Title", ""),
        record.get("Abstract", ""),
        record.get("Author Keywords", ""),
        record.get("Index Keywords", ""),
    ]).casefold()
    score = 0
    for needle, weight in KEYWORD_WEIGHTS.items():
        if needle in text:
            score += weight
    score += 4 * len(query_hits)
    cited_by = record.get("Cited by") or "0"
    try:
        score += min(10, int(cited_by) // 10)
    except ValueError:
        pass
    if "s3-symmetry-counts" in query_hits:
        score += 5
    if "s4-ddgp-counting-limits" in query_hits:
        score += 5
    if "s6-local-symmetry" in query_hits:
        score += 4
    return score


def compact_record(record_id: str, rows: list[dict[str, str]]) -> dict[str, object]:
    first = rows[0]
    query_hits = {row["_query"] for row in rows}
    tags = tags_for(first)
    score = score_record(first, query_hits)
    return {
        "id": record_id,
        "score": score,
        "query_hits": sorted(query_hits),
        "tags_hint": tags,
        "title": clean_space(first.get("Title")),
        "year": clean_space(first.get("Year")),
        "venue": clean_space(first.get("Source title")),
        "doi": clean_space(first.get("DOI")),
        "cited_by": clean_space(first.get("Cited by")),
        "document_type": clean_space(first.get("Document Type")),
        "authors": clean_space(first.get("Authors")),
        "author_keywords": clean_space(first.get("Author Keywords")),
        "index_keywords": clean_space(first.get("Index Keywords")),
        "abstract": truncate(first.get("Abstract", "")),
    }


def write_batch(batch_no: int, batch: list[dict[str, object]]) -> None:
    path = BATCH_DIR / f"scopus-batch-{batch_no:02d}.md"
    result_path = RESULT_DIR / f"scopus-batch-{batch_no:02d}-screening.md"
    lines = [
        f"# Scopus Screening Batch {batch_no:02d}",
        "",
        f"Records: {len(batch)}",
        f"Result file to write: `references/screening/agent-results/{result_path.name}`",
        "",
        "Classify each record as `core`, `background`, `maybe`, or `reject`.",
        "Use the tag vocabulary from the reference-screen skill. Give one sentence",
        "of rationale for every `core` and `maybe` record; terse rationale is enough",
        "for `background` and `reject` records.",
        "",
    ]
    for rec in batch:
        lines.extend([
            f"## {rec['id']}. {rec['title']}",
            "",
            f"- Year: {rec['year']}",
            f"- Venue: {rec['venue']}",
            f"- DOI: {rec['doi'] or 'n/a'}",
            f"- Cited by: {rec['cited_by'] or '0'}",
            f"- Document type: {rec['document_type'] or 'n/a'}",
            f"- Query hits: {', '.join(rec['query_hits'])}",
            f"- Tags hint: {', '.join(rec['tags_hint']) if rec['tags_hint'] else 'none'}",
            f"- Authors: {rec['authors']}",
            f"- Author keywords: {rec['author_keywords'] or 'n/a'}",
            f"- Index keywords: {rec['index_keywords'] or 'n/a'}",
            "",
            f"Abstract: {rec['abstract'] or 'n/a'}",
            "",
        ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    SCREEN_DIR.mkdir(parents=True, exist_ok=True)
    BATCH_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    raw_rows = 0
    for path in sorted(SEARCH_DIR.glob("scopus-*.csv.csv")):
        qid = query_id(path)
        with path.open(newline="", encoding="utf-8-sig", errors="replace") as f:
            for row in csv.DictReader(f):
                raw_rows += 1
                row["_query"] = qid
                doi = clean_space(row.get("DOI")).casefold()
                key = doi or normalize_title(row.get("Title"))
                grouped[key].append(row)

    records = [
        compact_record(f"R{i:03d}", rows)
        for i, rows in enumerate(grouped.values(), start=1)
    ]
    records.sort(
        key=lambda r: (
            -int(r["score"]),
            -len(r["query_hits"]),
            -(int(r["cited_by"]) if str(r["cited_by"]).isdigit() else 0),
            str(r["year"]),
            str(r["title"]).casefold(),
        )
    )
    for i, record in enumerate(records, start=1):
        record["id"] = f"R{i:03d}"

    compact_json = SCREEN_DIR / "scopus-compact-records.json"
    compact_json.write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")

    summary = {
        "raw_rows": raw_rows,
        "unique_records": len(records),
        "duplicate_rows": raw_rows - len(records),
        "batch_size": BATCH_SIZE,
        "batch_count": math.ceil(len(records) / BATCH_SIZE),
    }
    (SCREEN_DIR / "scopus-compact-summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    for old in BATCH_DIR.glob("scopus-batch-*.md"):
        old.unlink()
    for batch_no, start in enumerate(range(0, len(records), BATCH_SIZE), start=1):
        write_batch(batch_no, records[start:start + BATCH_SIZE])

    readme = [
        "# Scopus Screening Batches",
        "",
        f"Raw rows: {summary['raw_rows']}",
        f"Unique records: {summary['unique_records']}",
        f"Duplicate rows removed: {summary['duplicate_rows']}",
        f"Batch size: {summary['batch_size']}",
        f"Batch count: {summary['batch_count']}",
        "",
        "Generated from `references/searches/scopus-*.csv.csv`.",
        "Each batch should be screened independently into the matching file under",
        "`references/screening/agent-results/`.",
    ]
    (BATCH_DIR / "README.md").write_text("\n".join(readme), encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
