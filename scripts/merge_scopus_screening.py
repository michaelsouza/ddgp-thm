#!/usr/bin/env python3
"""Merge per-batch Scopus screening results into shortlist files."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCREEN_DIR = ROOT / "references" / "screening"
RESULT_DIR = SCREEN_DIR / "agent-results"
COMPACT_JSON = SCREEN_DIR / "scopus-compact-records.json"

CLASS_ORDER = {
    "core": 0,
    "maybe": 1,
    "background": 2,
    "reject": 3,
}


def clean(value: str | None) -> str:
    return re.sub(r"\s+", " ", value or "").strip()


def parse_table_line(line: str) -> tuple[str, dict[str, str]] | None:
    if not line.startswith("| R"):
        return None
    cells = [clean(cell) for cell in line.strip().strip("|").split("|")]
    if len(cells) < 8:
        return None
    rid, classification, doi, title, year, venue, tags, rationale = cells[:8]
    if not re.fullmatch(r"R\d{3}", rid):
        return None
    return rid, {
        "classification": classification.casefold(),
        "doi": doi,
        "title": title,
        "year": year,
        "venue": venue,
        "tags": tags,
        "rationale": rationale,
    }


def parse_section(lines: list[str], start: int) -> tuple[str, dict[str, str], int] | None:
    heading = lines[start]
    match = re.match(r"^##\s+(R\d{3})(?:\.\s*(.*))?$", heading)
    if not match:
        return None
    rid = match.group(1)
    heading_title = clean(match.group(2))
    fields: dict[str, str] = {}
    i = start + 1
    while i < len(lines) and not lines[i].startswith("## "):
        line = lines[i]
        bullet = re.match(r"^-\s*([^:]+):\s*(.*)$", line)
        if bullet:
            key = clean(bullet.group(1)).casefold().replace(" ", "_")
            fields[key] = clean(bullet.group(2))
        i += 1
    parsed = {
        "classification": fields.get("classification", "").casefold(),
        "doi": fields.get("doi", ""),
        "title": fields.get("title", heading_title),
        "year": fields.get("year", ""),
        "venue": fields.get("venue", ""),
        "tags": fields.get("tags", ""),
        "rationale": fields.get("rationale", ""),
    }
    return rid, parsed, i


def parse_result_file(path: Path) -> dict[str, dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    out: dict[str, dict[str, str]] = {}
    i = 0
    while i < len(lines):
        table = parse_table_line(lines[i])
        if table:
            rid, parsed = table
            out[rid] = parsed
            i += 1
            continue
        section = parse_section(lines, i)
        if section:
            rid, parsed, next_i = section
            out[rid] = parsed
            i = next_i
            continue
        i += 1
    return out


def load_screening() -> dict[str, dict[str, str]]:
    screening: dict[str, dict[str, str]] = {}
    for path in sorted(RESULT_DIR.glob("scopus-batch-*-screening.md")):
        for rid, parsed in parse_result_file(path).items():
            parsed["source_file"] = path.name
            screening[rid] = parsed
    return screening


def merge_records() -> list[dict[str, object]]:
    compact = json.loads(COMPACT_JSON.read_text(encoding="utf-8"))
    screening = load_screening()
    merged = []
    missing = []
    bad_class = []
    for record in compact:
        rid = record["id"]
        parsed = screening.get(rid)
        if not parsed:
            missing.append(rid)
            continue
        classification = parsed.get("classification", "").casefold()
        if classification not in CLASS_ORDER:
            bad_class.append((rid, classification))
        merged.append({
            **record,
            "classification": classification,
            "screening_tags": parsed.get("tags") or "; ".join(record.get("tags_hint", [])),
            "rationale": parsed.get("rationale", ""),
            "source_file": parsed.get("source_file", ""),
        })
    if missing or bad_class:
        details = []
        if missing:
            details.append(f"missing={missing}")
        if bad_class:
            details.append(f"bad_class={bad_class}")
        raise SystemExit("Invalid screening merge: " + "; ".join(details))
    return merged


def record_sort_key(record: dict[str, object]) -> tuple:
    cited = str(record.get("cited_by", "0"))
    cited_i = int(cited) if cited.isdigit() else 0
    return (
        CLASS_ORDER[str(record["classification"])],
        -int(record.get("score", 0)),
        -len(record.get("query_hits", [])),
        -cited_i,
        str(record.get("year", "")),
        str(record.get("title", "")).casefold(),
    )


def metadata_line(record: dict[str, object]) -> str:
    doi = record.get("doi") or "n/a"
    cited = record.get("cited_by") or "0"
    queries = ", ".join(record.get("query_hits", []))
    return (
        f"- DOI: {doi}\n"
        f"- Year: {record.get('year') or 'n/a'}\n"
        f"- Venue: {record.get('venue') or 'n/a'}\n"
        f"- Cited by: {cited}\n"
        f"- Query hits: {queries or 'n/a'}\n"
        f"- Tags: {record.get('screening_tags') or 'none'}\n"
        f"- Rationale: {record.get('rationale') or 'n/a'}"
    )


def write_shortlist(records: list[dict[str, object]]) -> None:
    counts = Counter(str(r["classification"]) for r in records)
    lines = [
        "# Scopus Screening Shortlist",
        "",
        "Generated from compacted Scopus exports and per-batch screening results.",
        "",
        "## Summary",
        "",
        f"- Total screened records: {len(records)}",
        f"- Core: {counts['core']}",
        f"- Maybe: {counts['maybe']}",
        f"- Background: {counts['background']}",
        f"- Reject: {counts['reject']}",
        "",
        "## Core",
        "",
        "Download or verify full text for all 25 core records. They are likely",
        "citation candidates for the camera-ready article.",
        "",
    ]
    for record in sorted([r for r in records if r["classification"] == "core"], key=record_sort_key):
        lines.extend([
            f"### {record['id']}. {record['title']}",
            "",
            metadata_line(record),
            "",
        ])
    lines.extend([
        "## Maybe",
        "",
        "Check full text only if the core set leaves a gap.",
        "",
    ])
    for record in sorted([r for r in records if r["classification"] == "maybe"], key=record_sort_key):
        lines.extend([
            f"### {record['id']}. {record['title']}",
            "",
            metadata_line(record),
            "",
        ])
    lines.extend([
        "## Background",
        "",
        "Useful context, but probably not needed in the short article unless a",
        "specific narrative gap appears.",
        "",
    ])
    for record in sorted([r for r in records if r["classification"] == "background"], key=record_sort_key):
        lines.extend([
            f"### {record['id']}. {record['title']}",
            "",
            metadata_line(record),
            "",
        ])
    (SCREEN_DIR / "shortlist.md").write_text("\n".join(lines), encoding="utf-8")


def write_rejected(records: list[dict[str, object]]) -> None:
    lines = [
        "# Scopus Screening Rejected Records",
        "",
        "Rejected records are mostly acronym collisions, generic application papers,",
        "or distance-geometry material outside the DDGP/DMDGP rank-count scope.",
        "",
    ]
    for record in sorted([r for r in records if r["classification"] == "reject"], key=record_sort_key):
        lines.extend([
            f"## {record['id']}. {record['title']}",
            "",
            metadata_line(record),
            "",
        ])
    (SCREEN_DIR / "rejected.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    records = merge_records()
    write_shortlist(records)
    write_rejected(records)
    counts = Counter(str(r["classification"]) for r in records)
    print(json.dumps({
        "total": len(records),
        "core": counts["core"],
        "maybe": counts["maybe"],
        "background": counts["background"],
        "reject": counts["reject"],
        "shortlist": str(SCREEN_DIR / "shortlist.md"),
        "rejected": str(SCREEN_DIR / "rejected.md"),
    }, indent=2))


if __name__ == "__main__":
    main()
