# DMDGP Symmetry/Counting Citation Report

## 1. Papers Reviewed

- `liberti2008branch`: original Branch-and-Prune presentation for the molecular distance geometry setting; supports discrete two-branch search and complete enumeration by BP, not the later symmetry/counting theorems.
- `liberti2011number`: COCOA conference version of the DMDGP power-of-two solution-count result; useful historically, but superseded by the 2014 journal version for theorem claims.
- `mucherino2012exploiting`: symmetry-driven BP (`symBP`); supports binary branch-code flipping, symmetry-set use, and reconstructing/enumerating solutions from symmetry information after one BP solution.
- `liberti2014number`: definitive primary theorem source for K-DMDGP partial-reflection structure, generic/almost-always power-of-two counts, and the warning that extending the result from DMDGP to nonconsecutive DDGP is nontrivial.
- `fidalgo2018symmetry`: symmetry-based splitting/covering strategy; cite only if the article discusses decomposing DMDGP instances into subinstances or symmetry covers.
- `lavor2021optimality`: optimal algorithm for finding DMDGP symmetry vertices; cite for efficient/a-priori symmetry-set detection, not as the primary source of the count theorem.
- `goncalves2021new`: SBBU algorithm; supports symmetry-based build-up from nested pruning-edge subproblems and gives a modern K-DMDGP formulation with partial reflections.

## 2. Citation-Supported Claims We Can Safely Make

- BP discretizes the relevant molecular/DMDGP search into a binary tree and can enumerate all feasible solutions: cite `liberti2008branch`; for modern K-DMDGP framing, cite `goncalves2021new` if needed.
- In DMDGP/K-DMDGP, consecutiveness of the predecessor order is essential for suffix/partial-reflection symmetry arguments: cite `liberti2014number`.
- For feasible K-DMDGP instances with exact real edge weights, the number of incongruent realizations is a power of two with probability one / outside a measure-zero exceptional set: cite `liberti2014number` as the primary source; avoid relying on `liberti2011number` for the camera-ready theorem statement.
- Valid DMDGP branchings induce partial reflections, and the solution set can be organized by binary branch codes and reflection operations: cite `liberti2014number`; add `mucherino2012exploiting` for the algorithmic branch-vector/flipping presentation.
- The exact-distance symmetry theory should not be stated for general nonconsecutive DDGP without qualification; `liberti2014number` gives a DDGP-not-DMDGP counterexample, and the current article's local/generic restriction is appropriate.
- DMDGP symmetry vertices/sets can be detected optimally in `O(|V| + |E_p|)`: cite `lavor2021optimality`.
- Symmetry-aware build-up can find a first realization by solving nested pruning-edge subproblems and using partial reflections: cite `goncalves2021new`.
- Symmetry-based splitting into subinstances is supported if the article later discusses decomposition or covers: cite `fidalgo2018symmetry`; otherwise omit it from the current text.

## 3. Existing Article Citations

- Intro, lines 47-53: keep `liberti2008branch` for BP, and keep `mucherino2012exploiting`/`goncalves2021new` for symmetry-aware algorithms. Add `liberti2014number` for the theorem-level symmetry/counting claim. Prefer changing the count sentence to cite `liberti2014number` directly, e.g. "Those symmetries underlie generic power-of-two solution counts \cite{liberti2014number}, symmetry-driven enumeration \cite{mucherino2012exploiting}, optimal symmetry detection \cite{lavor2021optimality}, and symmetry-based build-up algorithms \cite{goncalves2021new}."
- Intro, current use of `liberti2011number`: change to `liberti2014number` wherever the claim is theorem-level. Keep `liberti2011number` only if the article explicitly wants to acknowledge the conference predecessor.
- Lines 78-80: `mucherino2012exploiting` is appropriate for binary branch-code/symmetry-vector usage; no change required.
- Scope/genericity, line 575: replace `liberti2011number` with `liberti2014number` for the DMDGP almost-always counting qualification; keep `abud2024impossible` for the general-DDGP negative/counting caveat.
- Do not add `fidalgo2018symmetry` unless the article mentions symmetry splitting, subinstance covers, or decomposition.

## 4. BibTeX Metadata Issues Noticed

- `mucherino2012exploiting`: add DOI `10.1142/S0219720012420097`; local notes identify volume 10, number 3, article/page `1242009`.
- `liberti2014number`: add DOI `10.1016/j.dam.2013.01.020`.
- `fidalgo2018symmetry`: add DOI `10.1007/s10898-018-0610-9`.
- `lavor2021optimality`: add DOI `10.1007/s40314-021-01479-6`; fix author accent from `Andr{\^e}s` to `Andr{\'e}s`.
- `goncalves2021new`: add DOI `10.1007/s00453-021-00835-6`; protect title capitalization for `{K}` and `{DMDGP}`; consider encoding the author as `Gon{\c{c}}alves, Douglas S.`.
- `liberti2011number`: if retained, improve proceedings metadata: LNCS volume 6831, booktitle `Combinatorial Optimization, Constraints and Applications`, publisher Springer, pages 322--342.
- `liberti2008branch`: no DOI was visible in the local notes; do not invent one. The `publisher={Wiley Online Library}` field looks like a database label rather than a publisher.
