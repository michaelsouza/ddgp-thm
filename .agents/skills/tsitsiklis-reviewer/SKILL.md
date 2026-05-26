---
name: tsitsiklis-reviewer
description: Meticulous academic review of a user-specified LaTeX manuscript with mathematical content, based on John N. Tsitsiklis's writing recommendations. Use when reviewing an indicated article, paper folder, main .tex file, or LaTeX section bundle for structure, mathematical flow, notation, prose, typesetting, figures, and compilation-halting issues.
---

# Tsitsiklis Reviewer Skill

A systematic, publication-grade review workflow for LaTeX manuscripts with substantial mathematical content. The skill converts John N. Tsitsiklis's essay, *A Few Tips on Writing Papers with Mathematical Content*, into a manuscript-specific review protocol.

## Required Source

Before reviewing any manuscript, read the bundled reference [references/tsitsiklis2009few.md](references/tsitsiklis2009few.md), resolved relative to this skill folder. This file is the source of truth for the review standard. If the bundled file is missing, search the project for `tsitsiklis2009few.md`; if it cannot be found, tell the user and continue only with the summarized checklist below.

When using subagents, every subagent prompt must explicitly instruct the subagent to read the bundled `references/tsitsiklis2009few.md` before reading its assigned manuscript files.

## Input Scope

Review the article indicated by the user. The input may be a single `main.tex`, a paper folder, or a set of LaTeX section files. Do not assume `article/main.tex` unless the user names it or it is the only unambiguous manuscript entry point.

If the input is a main LaTeX file, inspect its `\input`, `\include`, bibliography, figure, table, and appendix structure to identify the review scope. If the input is a folder, find the likely main file with `rg --files -g '*.tex'` and by checking for `\documentclass`.

The default output is a review report, not direct edits. If the user asks to implement changes, apply only defensible fixes and preserve the author's mathematical intent.

## Review Standard

Use Tsitsiklis's full essay as the authoritative standard. The checklist below is a navigation aid, not a replacement for the essay.

### 1. Before-Writing Discipline

- Identify the intended audience implied by the paper and flag mismatches in assumed background.
- State the paper's main point in one or two sentences; if this is hard, the manuscript likely needs sharper framing.
- Check whether the table of contents, theorem order, notation, and terminology support that main point.
- Look for inconsistent spelling, hyphenation, terminology, graph-theoretic names, algorithm names, and symbol choices.

### 2. Document Structure

- The abstract should be crisp, direct, and light on notation. It should say what the paper considers and what it shows, without a long history or promotional language.
- The introduction should frame the problem early, then give motivation, context, literature, a preview of results, and an organization paragraph.
- Related work should compare with the closest work, not merely list citations.
- Sections and subsections should have one clear theme, useful headings, and short orientation text at the beginning.
- Long units should be split into digestible modules: definitions, assumptions, examples, lemmas, theorems, proofs, figures, and transitions.
- The conclusion should not merely repeat the abstract. It should emphasize ideas that became clear only after the technical development and should identify plausible next questions.
- Appendices may hold lengthy or routine proofs, but the main text must remain self-contained. Main-text notation and invoked results should not live only in an appendix.

### 3. Prose Style

- Prefer consistent first-person plural active voice: "We prove," "We define," "We use."
- Break long sentences, especially those carrying mathematical dependencies.
- Favor parallel constructions when presenting analogous cases, hypotheses, or conclusions.
- Use passive voice only deliberately, usually to omit an elementary argument that a competent reader can reconstruct; otherwise provide a hint or citation.
- Remove unnecessary words, but keep helpful parsing words such as "then" after an initial "if" clause and "that" after "assume" or "suppose" when clarity improves.
- Check pronouns such as "it," "this," "that," "which," and "they" for unambiguous antecedents.
- Avoid claims such as "easy," "trivial," "obvious," or "clearly" when they may irritate or slow a reader.
- Apply US-style punctuation and abbreviation conventions unless the manuscript or venue explicitly requires another style.

### 4. Mathematical Flow

- Every theorem, lemma, proposition, definition, and assumption should appear before the proof or discussion that depends on it.
- Formal statements should be short and crisp. Move terminology, notation, and consequence discussion outside the statement when possible.
- Define symbols close to first use, either immediately before use or immediately after use with a clear "where" clause.
- Do not introduce notation for expressions used only once or twice unless it genuinely improves readability.
- Introduce mnemonic labels for important assumptions, equations, and conditions instead of generic references.
- Proofs should flow forward: facts first, consequences next, conclusion last. Avoid proofs that alternate between backward planning and forward derivation.
- If a proof summary is useful, place it before the proof and distinguish motivation from rigorous argument.
- Avoid "magic trick" exposition: do not reveal a theorem only after proving it.
- Use examples, counterexamples, and figures to explain why assumptions are present or why they are stronger than necessary.
- Minimize reader effort by reminding the reader of distant definitions or equation meanings when needed.
- Flag ambiguous quantifiers and ambiguous asymptotic notation. State whether constants depend on parameters.

### 5. Mathematical Prose

- Sentences containing formulas should read as ordinary sentences.
- Do not start a sentence with a mathematical symbol.
- Avoid raw `\forall` and `\exists` in prose. Write "for every" and "there exists" instead.
- Do not say "the optimal solution," "the embedding," or similar unique language unless uniqueness has been proved.
- Punctuate displayed equations as part of the sentence.
- Ensure displayed equations have enough surrounding prose to tell the reader why they appear.

### 6. Low-Level LaTeX And Typesetting

- Use nonbreaking spaces before references and citations: `Theorem~\ref{...}`, `Eq.~\eqref{...}`, and `Author~\cite{...}`.
- Capitalize named references consistently: "Theorem," "Lemma," "Section," and "Equation" when referring to numbered items.
- Use proper delimiters: `\lvert ... \rvert`, `\lVert ... \rVert`, `\mid`, `\bigl`, `\bigr`, and related spacing commands where they improve readability.
- Avoid inline fractions that disrupt line spacing; prefer `(x+2)/(x+3)` when adequate.
- Set operators such as `\inf`, `\sup`, `\max`, `\argmin`, and `\rank` in roman/operator style.
- Add small spaces in integrals and dense formulas, such as `\,dx`.
- Consider macros for complicated repeated notation.
- Check LaTeX spacing after abbreviations such as `et~al.\`, `i.e.,`, `e.g.,`, and `cf.`
- Check display alignment, equation breaking, quote marks, bibliography punctuation, and figure/table captions.

### 7. Figures, Tables, Examples, And Counterexamples

- Captions should be self-contained enough for a reader to understand the figure without hunting through the text.
- Figures should illustrate results or intuition, not merely decorate.
- Examples should be concrete and checkable.
- Counterexamples should show why a hypothesis matters or explain when a hypothesis is stronger than necessary.
- Tables should have clear column meanings, units, and notation reminders.

## Review Workflow

1. Confirm the manuscript entry point or infer it only when unambiguous.
2. Read the bundled `references/tsitsiklis2009few.md`.
3. Map the manuscript structure: main file, included sections, appendices, bibliography, figures, and macros.
4. Compile the manuscript if a project-local compile command exists or if compilation can be done safely with installed tools. Report compile errors separately from writing issues.
5. Review high-level structure before line-level edits.
6. Review notation, theorem order, proof flow, examples, figures, and low-level LaTeX.
7. Distinguish objective problems from taste preferences. Do not recommend changes that would weaken mathematical precision.
8. Save the consolidated report in `reviews/review-YYYYMMDD-HHMMSS.md` using the local timestamp.

## Multi-Agent Protocol

Use subagents for long manuscripts, multi-file manuscripts, or when the user asks for a careful review. Split by manuscript structure, not by hard-coded filenames.

Every subagent must:

- read the bundled `references/tsitsiklis2009few.md` first;
- review only its assigned files plus any shared notation/macros needed to understand them;
- cite the relevant Tsitsiklis section or principle in each issue;
- provide exact file/line references and actionable replacements;
- avoid editing files directly unless the user asked for implementation rather than review.

Suggested roles:

1. **Front Matter And Framing Reviewer**
   - Scope: abstract, introduction, related work, contributions, organization, conclusion.
   - Focus: immediate framing, audience, main point, literature comparison, conclusion value, active voice.

2. **Model And Notation Reviewer**
   - Scope: definitions, assumptions, notation sections, problem formulation, preliminaries.
   - Focus: symbol declarations, terminology consistency, theorem prerequisites, self-contained main text, unnecessary notation.

3. **Theory And Proof Reviewer**
   - Scope: lemmas, propositions, theorems, proofs, appendices containing proof material.
   - Focus: short statements, forward proof flow, dependency order, quantifiers, asymptotic dependencies, consequence placement.

4. **Examples, Figures, Tables, And Typesetting Reviewer**
   - Scope: examples, experiments if present, figures, tables, captions, equations, citations, bibliography style.
   - Focus: concrete checkability, caption self-containment, assumption-motivating examples, display punctuation, delimiter and spacing quality.

For small manuscripts, one agent may perform all roles sequentially.

## Issue Format

Every issue must be specific and actionable:

```markdown
### Issue N: Short Descriptive Title
* **Severity**: Critical | Major | Minor | Suggestion
* **File & Line**: `path/to/file.tex` (Line N)
* **Tsitsiklis Basis**: Section or principle from the bundled `references/tsitsiklis2009few.md`
* **Problem**: Explain what slows, misleads, or burdens the reader.
* **Original LaTeX**:
  ```latex
  ...
  ```
* **Recommended Replacement**:
  ```latex
  ...
  ```
* **Rationale**: Explain why the replacement is better and whether it changes meaning.
```

Use "Critical" for compile failures, undefined references, broken math environments, or errors that can invalidate the manuscript. Use "Major" for structural, mathematical, or proof-flow problems. Use "Minor" for localized prose, notation, and LaTeX issues. Use "Suggestion" for optional improvements.

## Report Format

The consolidated report must use this structure:

```markdown
# Tsitsiklis-Style Review

* **Manuscript**: path or file set reviewed
* **Timestamp**: YYYY-MM-DD HH:MM:SS local time
* **Reference Standard**: bundled `references/tsitsiklis2009few.md`

## Executive Summary
[One short paragraph.]

## Critical Issues
[Compile failures or mathematically blocking problems first. Write "None found" if none.]

## Major Issues
[Structure, proof flow, notation, self-containedness.]

## Minor Issues
[Localized prose and typesetting.]

## Optional Suggestions
[Taste-level improvements.]

## Checks Performed
[Files read, compile command if any, scans used.]
```

Save the final report as `reviews/review-YYYYMMDD-HHMMSS.md`. Present the path and a concise summary to the user.
