# Handoff: DDGP Camera-Ready Rank-Count Article

## Purpose

This handoff is for a fresh agent/co-author to continue preparing the DDGP
rank-count theory as a camera-ready LaTeX article. The immediate next task is
to strengthen the proof of the Feasible Branch Shift Theorem.

## Primary Tracking Artifacts

- Parent PRD issue: https://github.com/michaelsouza/ddgp-thm/issues/1
- Next recommended issue: https://github.com/michaelsouza/ddgp-thm/issues/2
- Child issue list:
  - #2 Harden the feasible branch shift proof
  - #3 Polish the camera-ready article narrative
  - #4 Verify the fully worked local example
  - #5 Normalize article-facing terminology and glossary
  - #6 Plan HITL literature searches for related work
  - #7 Screen exported references for article relevance
  - #8 Integrate selected citations into the LaTeX article
  - #9 Run final build and reviewer sanity pass

Do not duplicate the PRD or issue bodies; use the GitHub issues as the source of
scope and acceptance criteria.

## Current Local Artifacts

- Article source: `article/main.tex`
- Generated PDF: `article/main.pdf`
- Bibliography: `references/references.bib`
- Short source report: `research/ddgp-rank-count-theory-report-short-formal.md`
- Project-local reference workflow skills:
  - `.agents/skills/literature-search-planner/SKILL.md`
  - `.agents/skills/reference-screen/SKILL.md`
  - `.agents/skills/citation-integrator/SKILL.md`

The article currently builds cleanly with LaTeX/BibTeX and contains no TODOs or
unresolved citations as of this handoff.

Build command sequence:

```bash
cd article
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Key Decisions Already Made

- The article is in LaTeX at `article/main.tex`.
- The bibliography is canonical at `references/references.bib`.
- BibTeX keys should follow Google-Scholar-style naming:
  `firstauthorYYYYfirstmeaningfultitleword`.
- Use US spelling: `labeled`, not `labelled`, in article-facing prose.
- Canonical terms:
  - `cone generator`
  - `base generator`
  - `labeled violation`
  - `feasible branch shift space`
  - `generic mirror separation`
- Avoid old canonical terms:
  - `dependency-cone generator`
  - `base-clique closure generator`
  - `stabilizer space` as the main name for `\mathcal K_F`
- Large numerical experiments are internal validation, not a central paper
  result. The article should use only a small fully worked example for the
  reader/reviewer.
- Keep the proof in the main body while it remains compact. Use an appendix
  only if the proof expands into technical obstruction or bracket-separation
  machinery.
- Do not migrate Python identifiers or CSV schema fields from `labelled_*` to
  `labeled_*` unless a later issue explicitly asks for a code/data migration.

## Current Article Shape

`article/main.tex` contains:

- Introduction
- Local DDGP subproblems
- Graph-derived generators
- Labeled violations
- Feasible branch shifts
- Main results
- Proof of the Feasible Branch Shift Theorem
- A small worked example
- Scope and genericity
- Conclusion

The main theorem statement is:

```tex
\Xi_F=s^\ast\oplus \mathcal K_F=s^\ast\oplus M(\ker V).
```

The rank-count formula is:

```tex
|\Xi_F|
=
2^{
\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}
-
\operatorname{rank}(V)
}.
```

## Recommended Next Step

Start with issue #2. The current proof is compact but still has a prose-heavy
step:

- zero-labeled presentations preserve active edges by grouping selected
  generators by mirror clique.

The best improvement is to add a short lemma before the theorem or proof,
roughly:

```tex
V\alpha=0 \implies M\alpha \text{ preserves all active edges.}
```

Then revise the theorem proof to invoke that lemma directly. Keep the
rank-nullity proof simple and central. Do not reintroduce the obstruction
quotient unless it becomes necessary.

## Suggested Skills

- `diagnose`: if LaTeX build, BibTeX, or theorem/example consistency checks fail.
- `literature-search-planner`: when starting issue #6; it is HITL because the
  user must log into Scopus/Web of Science and export CSV files.
- `reference-screen`: after CSV exports are present for issue #7.
- `citation-integrator`: after selected papers are downloaded/converted for
  issue #8.
- `grill-with-docs`: if terminology or proof structure needs to be challenged
  against the project glossary.
- `to-issues`: only if the parent PRD needs further breakdown.

## Cautions

- The working tree is intentionally dirty with article creation, terminology
  renames, generated skills, and article-facing spelling changes. Do not revert
  unrelated changes.
- Some old `labelled_*` names may remain in scripts or CSV files by design.
  Do not treat those as article-facing spelling bugs.
- Before committing, inspect generated LaTeX auxiliaries under `article/` and
  decide which generated files should be tracked. The source and PDF currently
  exist locally; auxiliary files may also exist from builds.
