# ddgp-thm

Research notes and glossary for the study of the Distance Geometry Problem (DGP) and related topics in mathematics and combinatorial optimization.

## Structure

```
concepts/       — glossary entries (one file per concept)
references/     — source articles (.pdf + .md notes)
research/       — evolving hypothesis journals, proofs, and the unified rank-count theory report
scripts/        — tooling
```

## Article review PDF

To build a review copy of the LaTeX article with physical source-line numbers
overlaid on the rendered mathematical PDF, run:

```bash
scripts/compile_numbered_article.sh
```

This generates `article/main_numbered.pdf`. The red margin labels are physical
line numbers from each `article/sections/*.tex` file, so numbering restarts for
each file included through `\input{sections/...}`. The generated numbered PDF and
intermediate files are ignored by git.

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

## Experimental tooling

The current DDGP counting experiments use:

- `scripts/ddgp_generate_single_pruning.py` to generate one-pruning-edge instances.
- `scripts/ddgp_generate_multi_pruning.py` to generate multi-pruning-edge instances.
- `scripts/ddgp_count_solutions.py` to enumerate exact counts for validation.
- `scripts/ddgp_rank_count.py` to predict local base and valid counts by the active-edge labeled-violation rank formula without enumerating branch strings.

The geometric generator and exact enumerator currently support `K=2` and `K=3`.

The current proof formulation for the rank count is in
`research/labeled-violation-rank-proof.md`. For `K >= 2`, the generic theorem
is formulated through the labeled presentation property and the labeled
obstruction quotient. A unified, article-style presentation of the complete theory,
proof architecture, F_2 rank-count algorithm, and experimental verification is available
in `research/ddgp-rank-count-theory-report.md`.

## Decisions

| Decision | Choice | Reason | Date |
|---|---|---|---|
| Concept file naming | English slugs | Consistent with technical terminology in source articles | 2026-05-20 |
| Concept scope | Broad (math + CS) | Articles span multiple areas; no premature narrowing | 2026-05-20 |
| Index management | Auto-generated from frontmatter | Avoids manual sync; single source of truth is each concept file | 2026-05-20 |
| Index grouping | First tag | Simple heuristic; author controls grouping by ordering tags | 2026-05-20 |
| Frontmatter parser | `python-frontmatter` | Robust YAML parsing; avoids hand-rolled regex | 2026-05-20 |
| Concept structure | Intuition + Formal definition | Two audiences: accessible reading and mathematical precision | 2026-05-20 |
| Solution-counting seed | Fix first `K+1` vertices | Removes the trivial global reflection; binary decisions start at vertex `K+2` | 2026-05-20 |
