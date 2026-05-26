# DDGP Camera-Ready Literature Search Plan

Date: 2026-05-25

## Search Objective

Find references for a camera-ready article on a generic rank-count theory for
local DDGP solution sets. The search should support four citation needs:

1. Position the work in distance geometry, DDGP, and DMDGP.
2. Cite Branch-and-Prune, discretization, and lateration background.
3. Cite DMDGP partial-reflection symmetry and solution-counting results.
4. Cite DDGP-specific limitations showing why a new local/generic counting
   theory is needed.

This is a human-in-the-loop search. Run the queries in Scopus and Web of
Science, export CSV records, and let the `reference-screen` skill rank the
results from titles, abstracts, keywords, venues, years, DOI, and citation data.

## Seed Papers

Use these as starting points for query refinement and citation chaining:

- Lavor, Liberti, Maculan, and Mucherino, "The discretizable distance geometry
  problem", 2012.
- Liberti, Masson, Lee, Lavor, and Mucherino, "On the number of solutions of
  the discretizable molecular distance geometry problem".
- Gonçalves, Lavor, Liberti, and Souza, "A new algorithm for the K-DMDGP
  subclass of distance geometry problems with exact distances", 2021.
- Gonçalves and Mucherino, "A new symmetry-based build-up algorithm for the
  discretizable molecular distance geometry problem", 2021.
- Abud, Alencar, Lavor, Liberti, and Mucherino, "An impossible combinatorial
  counting method in distance geometry", 2024.

## Search Facets

Problem terms:

- distance geometry
- discretizable distance geometry problem
- DDGP
- discretizable molecular distance geometry problem
- DMDGP
- molecular distance geometry

Algorithmic terms:

- Branch-and-Prune
- branch and prune
- lateration
- K-lateration
- sphere intersection
- pruning edge
- discretization edge

Symmetry/counting terms:

- partial reflection
- reflection symmetry
- symmetry group
- solution count
- number of solutions
- power of two
- combinatorial counting
- fixed-parameter tractability

Secondary mathematical terms, for exploratory searches only:

- genericity
- algebraic exceptional set
- rank over GF(2)
- binary branch code
- graph generator

## Scopus Queries

Run broad queries first. If a query returns more than roughly 300 records,
restrict by subject area to Mathematics, Computer Science, Engineering, or
Biochemistry/Genetics/Molecular Biology, but do not restrict by year initially.

### S1 Core DDGP/DMDGP

```text
TITLE-ABS-KEY(
  "discretizable distance geometry problem" OR
  "discretizable molecular distance geometry problem" OR
  DDGP OR DMDGP
)
```

### S2 Branch-and-Prune Background

```text
TITLE-ABS-KEY(
  ("distance geometry" OR "molecular distance geometry") AND
  ("branch-and-prune" OR "branch and prune" OR "Branch-and-Prune")
)
```

### S3 DMDGP Symmetry and Solution Counts

```text
TITLE-ABS-KEY(
  (DMDGP OR "discretizable molecular distance geometry problem") AND
  (symmetry OR symmetries OR "partial reflection" OR "partial reflections" OR
   "solution count" OR "number of solutions" OR "power of two")
)
```

### S4 DDGP Counting and Limitations

```text
TITLE-ABS-KEY(
  (DDGP OR "discretizable distance geometry problem") AND
  ("solution count" OR "number of solutions" OR "combinatorial counting" OR
   impossibility OR "fixed-parameter tractability" OR FPT)
)
```

### S5 Discretization and Lateration

```text
TITLE-ABS-KEY(
  ("distance geometry" OR "molecular distance geometry") AND
  (lateration OR "K-lateration" OR trilateration OR "sphere intersection") AND
  (discretization OR discretizable OR "Branch-and-Prune" OR "branch and prune")
)
```

### S6 Narrow Local Symmetry Query

```text
TITLE-ABS-KEY(
  ("distance geometry" OR DDGP OR DMDGP) AND
  ("partial reflection" OR "reflection symmetry" OR "symmetry-based") AND
  ("pruning edge" OR "pruning edges" OR "branch code" OR "binary tree")
)
```

### S7 Exploratory Rank/Code Query

Use this only after the core searches, because it may retrieve unrelated coding
or rigidity papers.

```text
TITLE-ABS-KEY(
  ("distance geometry" OR DDGP OR DMDGP) AND
  (rank OR "GF(2)" OR "finite field" OR "binary code" OR "branch code") AND
  (symmetry OR counting OR "solution count")
)
```

## Web of Science Queries

Run the same conceptual searches using `TS=`. Use `NEAR/x` where useful.

### W1 Core DDGP/DMDGP

```text
TS=(
  "discretizable distance geometry problem" OR
  "discretizable molecular distance geometry problem" OR
  DDGP OR DMDGP
)
```

### W2 Branch-and-Prune Background

```text
TS=(
  ("distance geometry" OR "molecular distance geometry") AND
  ("branch-and-prune" OR "branch and prune" OR "Branch-and-Prune")
)
```

### W3 DMDGP Symmetry and Solution Counts

```text
TS=(
  (DMDGP OR "discretizable molecular distance geometry problem") AND
  (symmetry OR symmetries OR "partial reflection" OR "partial reflections" OR
   "solution count" OR "number of solutions" OR "power of two")
)
```

### W4 DDGP Counting and Limitations

```text
TS=(
  (DDGP OR "discretizable distance geometry problem") AND
  ("solution count" OR "number of solutions" OR "combinatorial counting" OR
   impossibility OR "fixed-parameter tractability" OR FPT)
)
```

### W5 Discretization and Lateration

```text
TS=(
  ("distance geometry" OR "molecular distance geometry") AND
  (lateration OR "K-lateration" OR trilateration OR "sphere intersection") AND
  (discretization OR discretizable OR "Branch-and-Prune" OR "branch and prune")
)
```

### W6 Narrow Local Symmetry Query

```text
TS=(
  ("distance geometry" OR DDGP OR DMDGP) AND
  ("partial reflection" OR "reflection symmetry" OR "symmetry-based") AND
  ("pruning edge" OR "pruning edges" OR "branch code" OR "binary tree")
)
```

### W7 Citation-Chaining Query

Use Web of Science's "Cited References" or "Citation Network" tools on the seed
papers. Export papers that cite at least one of:

- "The discretizable distance geometry problem"
- "On the number of solutions of the discretizable molecular distance geometry
  problem"
- "An impossible combinatorial counting method in distance geometry"

## Export Instructions

For every query, export CSV or Excel with at least:

- title
- abstract
- authors
- year
- source title / venue
- DOI
- author keywords
- index keywords
- references, if available
- cited-by count, if available
- document type
- database source
- query identifier, if the database allows custom tags

Recommended filenames:

- `references/searches/scopus-s1-core-ddgp-2026-05-25.csv`
- `references/searches/scopus-s2-branch-prune-2026-05-25.csv`
- `references/searches/scopus-s3-symmetry-counts-2026-05-25.csv`
- `references/searches/scopus-s4-ddgp-counting-limits-2026-05-25.csv`
- `references/searches/scopus-s5-lateration-2026-05-25.csv`
- `references/searches/scopus-s6-local-symmetry-2026-05-25.csv`
- `references/searches/scopus-s7-rank-code-2026-05-25.csv`
- `references/searches/wos-w1-core-ddgp-2026-05-25.csv`
- `references/searches/wos-w2-branch-prune-2026-05-25.csv`
- `references/searches/wos-w3-symmetry-counts-2026-05-25.csv`
- `references/searches/wos-w4-ddgp-counting-limits-2026-05-25.csv`
- `references/searches/wos-w5-lateration-2026-05-25.csv`
- `references/searches/wos-w6-local-symmetry-2026-05-25.csv`
- `references/searches/wos-w7-citation-chain-2026-05-25.csv`

If a query returns many duplicate or irrelevant records, export anyway unless it
is larger than roughly 1000 records. The screening step can deduplicate by DOI
and normalized title.

## Inclusion Criteria for Screening

Keep papers that satisfy at least one criterion:

- Defines or surveys DGP, DDGP, DMDGP, or molecular distance geometry.
- Introduces or analyzes Branch-and-Prune, discretization orders, pruning edges,
  lateration, or sphere-intersection search.
- Studies DMDGP symmetries, partial reflections, symmetry-based build-up, or
  solution counts.
- Explains why DDGP differs from DMDGP, especially non-contiguous predecessor
  sets or order-recognition complexity.
- Discusses combinatorial versus metric/generic solution counting in distance
  geometry.
- Provides standard background on genericity or exceptional algebraic cases used
  in DDGP/DMDGP counting arguments.

## Exclusion Criteria for Screening

Exclude papers when the title/abstract clearly indicates:

- Only continuous global optimization for DGP, with no discretization,
  Branch-and-Prune, lateration, or symmetry connection.
- Sensor-network localization or EDM completion without a useful DDGP/DMDGP,
  lateration, or symmetry angle.
- Protein modeling papers that use distance geometry only as a black-box
  application and do not discuss the mathematical/combinatorial structure.
- Generic rigidity theory papers with no link to discretizable orders,
  lateration, solution counting, or reflection symmetries.
- Papers not in English unless they are foundational and metadata is sufficient.

## Screening Priorities

For `reference-screen`, rank records into:

- `must-cite`: foundational DDGP/DMDGP, BP, symmetry, solution-counting, or DDGP
  impossibility papers.
- `likely-cite`: directly related algorithms, surveys, or theoretical papers.
- `background-only`: useful context but probably not cited in the short article.
- `exclude`: unrelated after title/abstract screening.

Prefer a compact final bibliography. The article likely needs a focused set of
roughly 8-15 references, not an extensive survey.
