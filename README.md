# ddgp-thm

Research notes and glossary for the study of the Distance Geometry Problem (DGP) and related topics in mathematics and combinatorial optimization.

## Structure

```
concepts/       — glossary entries (one file per concept)
references/     — source articles (.pdf + .md notes)
scripts/        — tooling
```

## Glossary conventions

### Creating a concept

Copy the template and name the file using an English slug:

```bash
cp concepts/_template.md concepts/euclidean-distance.md
```

### File format

Each concept file uses YAML frontmatter followed by two fixed sections:

```markdown
---
term: Euclidean Distance       # display name (may include spaces, accents)
tags: [metric, geometry]       # first tag is used for grouping in the index
sources: [goncalves2021new]    # cite keys matching files in references/
---

## Intuition

Free-text explanation: what the concept means, why it matters,
and how it connects to related ideas.

## Formal definition

$$
d(x, y) = \sqrt{\sum_{i=1}^{n} (x_i - y_i)^2}
$$

Where $x, y \in \mathbb{R}^n$.
```

### Regenerating the index

```bash
source .venv/bin/activate
python3 scripts/build-index.py
```

`concepts/index.md` is auto-generated — do not edit it directly.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Decisions

| Decision | Choice | Reason | Date |
|---|---|---|---|
| Concept file naming | English slugs | Consistent with technical terminology in source articles | 2026-05-20 |
| Concept scope | Broad (math + CS) | Articles span multiple areas; no premature narrowing | 2026-05-20 |
| Index management | Auto-generated from frontmatter | Avoids manual sync; single source of truth is each concept file | 2026-05-20 |
| Index grouping | First tag | Simple heuristic; author controls grouping by ordering tags | 2026-05-20 |
| Frontmatter parser | `python-frontmatter` | Robust YAML parsing; avoids hand-rolled regex | 2026-05-20 |
| Concept structure | Intuition + Formal definition | Two audiences: accessible reading and mathematical precision | 2026-05-20 |
