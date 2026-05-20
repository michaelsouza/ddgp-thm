#!/usr/bin/env python3
"""Regenerate concepts/index.md from frontmatter of all concept files."""

import frontmatter
from pathlib import Path

CONCEPTS_DIR = Path(__file__).parent.parent / "concepts"
INDEX_FILE = CONCEPTS_DIR / "index.md"


def collect_entries() -> list[dict]:
    entries = []
    for path in sorted(CONCEPTS_DIR.glob("*.md")):
        if path.name.startswith("_") or path.name == "index.md":
            continue
        post = frontmatter.load(path)
        if not post.get("term"):
            continue
        entries.append(
            {
                "term": post["term"],
                "tags": post.get("tags") or [],
                "sources": post.get("sources") or [],
                "file": path.name,
            }
        )
    return entries


def build_index(entries: list[dict]) -> str:
    lines = [
        "# Concepts index\n",
        f"_Auto-generated. {len(entries)} concept(s). Edit individual files, not this index._\n",
    ]

    # Group by first tag (or "uncategorized")
    groups: dict[str, list[dict]] = {}
    for entry in entries:
        key = entry["tags"][0] if entry["tags"] else "uncategorized"
        groups.setdefault(key, []).append(entry)

    for group in sorted(groups):
        lines.append(f"\n## {group}\n")
        for entry in sorted(groups[group], key=lambda e: e["term"].lower()):
            tag_str = ", ".join(f"`{t}`" for t in entry["tags"]) if entry["tags"] else ""
            src_str = ", ".join(entry["sources"]) if entry["sources"] else ""
            detail = " — ".join(filter(None, [tag_str, src_str]))
            lines.append(f"- [{entry['term']}]({entry['file']})" + (f" ({detail})" if detail else ""))

    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    entries = collect_entries()
    INDEX_FILE.write_text(build_index(entries))
    print(f"index.md updated: {len(entries)} concept(s) indexed.")
