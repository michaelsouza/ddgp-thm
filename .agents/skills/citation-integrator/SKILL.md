---
name: citation-integrator
description: Integrate selected references into the DDGP camera-ready LaTeX article. Use when Codex needs to read downloaded paper Markdown/PDF notes, add or update BibTeX entries in references/references.bib, and place precise citations in the LaTeX article.
---

# Citation Integrator

Use this skill after selected papers have been downloaded and converted to
Markdown under `references/`. This is the only reference workflow skill that
should edit the LaTeX article.

## Required Paths

- Article source: `article/main.tex` by default, unless the user selects a
  different `.tex` file.
- Bibliography: `references/references.bib`.
- Paper notes: Markdown/PDF files under `references/`.

Create `references/references.bib` only when adding the first BibTeX entry.

## Workflow

1. Read the article section that needs citations.
2. Read only the selected reference Markdown files needed for that section.
3. Extract citation-relevant claims:
   - definitions and problem classes;
   - algorithmic baselines;
   - known symmetry/counting results;
   - limitations of prior DDGP/DMDGP theory.
4. Add or update BibTeX entries in `references/references.bib`.
5. Insert citations in the LaTeX source using existing citation style.
6. Keep claims narrow: cite only what the paper actually supports.
7. Prefer primary sources over survey or secondary references for theorem
   statements.

## Guardrails

- Use Google-Scholar-style BibTeX keys:
  `normalizedFirstAuthorSurnameYYYYfirstMeaningfulTitleWord`, for example
  `lavor2012discretizable`. Normalize accents to ASCII in keys.
- Do not cite a paper based only on title or abstract.
- Do not fabricate DOI, pages, issue, publisher, or BibTeX fields.
- Do not quote long passages from copyrighted papers.
- If metadata is incomplete, add the safest minimal BibTeX entry and mark the
  missing field in a nearby note or final response.
- Keep `references/references.bib` as the canonical bibliography file.
