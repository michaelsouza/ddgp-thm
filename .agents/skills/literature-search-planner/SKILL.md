---
name: literature-search-planner
description: Plan literature searches for the DDGP camera-ready article. Use when Codex needs to suggest Scopus, Web of Science, or scholarly web queries; define inclusion/exclusion criteria; map seed papers to search strings; or prepare instructions for exporting reference CSV files for later screening.
---

# Literature Search Planner

Use this skill to prepare database searches. Do not screen papers deeply and do
not edit the article.

## Workflow

1. Read the article draft or topic statement, plus any seed references in
   `references/`.
2. Identify search facets:
   - problem family: distance geometry, DDGP, DMDGP, discretizable distance geometry;
   - methods: Branch-and-Prune, partial reflections, symmetry, solution counting;
   - mathematics: rank over GF(2), graph generators, reflection groups, genericity.
3. Produce database-ready queries for Scopus and Web of Science.
4. Include broad and narrow query variants.
5. Specify export fields to request:
   - title;
   - abstract;
   - authors;
   - year;
   - venue;
   - DOI;
   - keywords;
   - references/cited-by count when available.
6. Define inclusion/exclusion criteria for the later screening skill.

## Output

Write concise Markdown suitable for saving under `references/searches/`, with:

- search objective;
- seed papers;
- Scopus queries;
- Web of Science queries;
- export instructions;
- inclusion/exclusion criteria.

Do not invent citation metadata. If exact metadata is needed, ask the user to
export it from Scopus/Web of Science or provide the downloaded CSV.
