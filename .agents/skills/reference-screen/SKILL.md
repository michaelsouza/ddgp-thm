---
name: reference-screen
description: Screen exported Scopus or Web of Science reference CSV files for the DDGP camera-ready article. Use when Codex needs to rank or filter references by title, abstract, keywords, DOI, venue, and year before the user downloads selected papers.
---

# Reference Screen

Use this skill after the user provides exported reference CSV files. Work from
titles, abstracts, keywords, venues, years, and DOI fields only unless full
papers are already present.

## Workflow

1. Locate user-provided CSV files, normally under `references/` or
   `references/searches/`.
2. Inspect headers before parsing; support Scopus and Web of Science naming
   differences.
3. Deduplicate records by DOI when available, otherwise normalized title.
4. Classify each record:
   - core: likely cite in the article;
   - background: useful context but not central;
   - maybe: needs full-text check;
   - reject: outside scope.
5. Tag relevance dimensions:
   - DDGP/DMDGP definitions;
   - Branch-and-Prune;
   - partial reflection symmetries;
   - solution counting;
   - distance geometry applications;
   - algebraic/genericity methods.
6. Explain each core and maybe decision in one sentence.

## Output

Create or update:

- `references/screening/shortlist.md`;
- `references/screening/rejected.md` when useful.

Keep the shortlist actionable: include DOI, title, year, venue, relevance tags,
and why the paper should or should not be downloaded.

Do not add BibTeX entries unless reliable citation metadata is already present.
Do not cite papers in the article; that belongs to `citation-integrator`.
