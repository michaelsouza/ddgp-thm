# DDGP Definition/Order/Core-Structure Citation Report

## 1. Papers reviewed.

- Mucherino, Lavor, Liberti, "The discretizable distance geometry problem" (`mucherino2012discretizable`).
- Cassioli, Gunluk, Lavor, Liberti, "Discretization vertex orders in distance geometry" (`cassioli2015discretization`).
- Abud, Alencar, Lavor, Liberti, Mucherino, "The K-discretization and K-incident graphs for discretizable Distance Geometry" (`abud2018k`).
- Abud, Alencar, Lavor, Liberti, Mucherino, "An impossible combinatorial counting method in distance geometry" (`abud2024impossible`).

## 2. Citation-supported claims we can safely make, with suggested BibTeX keys.

- The DGP can be stated as realizing a weighted graph in `R^K` with Euclidean distances matching edge weights. Suggested keys: `abud2024impossible`, `cassioli2015discretization`, `abud2018k`.
- The original DDGP paper defines a three-dimensional DDGP variant where each non-initial vertex has three adjacent predecessors, not necessarily the three immediate predecessors. Suggested key: `mucherino2012discretizable`.
- The K-dimensional DDGP order can be stated as an initial clique/order seed plus, for each later vertex, at least `K` adjacent predecessors in the order. Definitions vary in whether the selected predecessor set itself must induce a clique/positive simplex, so quote the exact convention being used. Suggested keys: `cassioli2015discretization`, `abud2024impossible`, `abud2018k`.
- DMDGP is the stricter subclass where the `K` adjacent predecessors are immediate and contiguous in the order; this contiguity implies strong clique/chain structure. Suggested keys: `cassioli2015discretization`, `abud2024impossible`.
- DDGP is broader than DMDGP because the predecessor set need not be contiguous; Mucherino et al. give explicit DDGP3 examples that are not DMDGP under any order. Suggested key: `mucherino2012discretizable`; for K-order language also cite `cassioli2015discretization`.
- BP/K-lateration uses the chosen predecessor set to generate generically at most two positions for the next vertex, then pruning edges/extra distances eliminate infeasible positions. Suggested keys: `mucherino2012discretizable`, `abud2024impossible`.
- Suitable DDGP orders are easier to recognize for fixed `K` than DMDGP/contiguous orders: Cassioli et al. prove CTOP is NP-complete for fixed `K`, while the TOP/DDGP-order difficulty is tied to the initial clique and is polynomial/FPT once `K` is fixed. Suggested key: `cassioli2015discretization`; `abud2024impossible` is useful for the later summary of this contrast.
- The order-existence hierarchy `CTOP subsetneq REOP subsetneq TOP` supports saying contiguous orders are more restrictive than repetition orders, which are more restrictive than DDGP orders. Suggested key: `cassioli2015discretization`.
- The K-discretization graph is order-independent, while the K-incident graph depends on a chosen discretization order and valid reference-vertex list; this is the closest reviewed prior art for graph structure around DDGP reference sets. Suggested key: `abud2018k`.
- The partition into discretization edges and pruning edges is explicit in Abud 2024, and Abud 2018 discusses how choices of reference vertices affect incident graphs/pruning structure. Suggested keys: `abud2024impossible`, `abud2018k`.
- DMDGP has graph-combinatorial solution counting through partial-reflection symmetry, but the same kind of full-instance, graph-only counting cannot extend to general DDGP because DDGP solution counts may depend on edge weights. Suggested key: `abud2024impossible`.
- Abud 2024 does not rule out special subclasses or local methods: it states that all `U_j` inducing cliques gives a combinatorial-DDGP sufficient condition without pruning edges, and it explicitly leaves room for methods using extra structure. Suggested key: `abud2024impossible`.
- The article's dependency cones, cone/base generators, labeled violation matrix, and rank-count theorem are not directly established in these reviewed papers. Cite `abud2018k`/`abud2024impossible` only for predecessor/reference/pruning context, and present the generator/violation machinery as this article's contribution.

## 3. Existing article citations that should remain/change/remove.

- Keep `mucherino2012discretizable` in the introduction for DDGP origin/BP generalization, but add `abud2024impossible` or `abud2018k` if the sentence is meant as the modern K-dimensional DDGP definition. Mucherino 2012 is DDGP3-focused.
- Be careful with the article's fixed `V_0 = {v_1, ..., v_{K+1}}` convention. Mucherino 2012 supports a four-vertex initial set in `R^3`; Cassioli 2015 and Abud 2024 state their order definitions with an initial `K`-clique. If the article keeps `K+1` fixed seeds to remove the initial reflection, present that as the article's convention rather than the common DDGP definition.
- Keep `cassioli2015discretization,abud2024impossible` for the sentence distinguishing general DDGP orders from DMDGP/contiguous orders.
- Keep `abud2024impossible` for the negative claim that full general DDGP solution counts cannot be graph-only/combinatorial.
- Add `abud2018k` near the local predecessor/reference-set setup or before the graph-derived generator section, because it is the reviewed paper most directly about K-discretization/K-incident structure.
- In the conclusion, qualify "combinatorial counting method for general DDGP instances" as "for fixed local active DDGP subproblems under the stated generic mirror-separation assumptions." Without that qualifier, it reads too close to the general counting method ruled out by `abud2024impossible`.
- No reviewed citation should be removed outright. Do not use `abud2024impossible` as support for a positive global DDGP counting theorem; use it as boundary/negative-result context.

## 4. Any missing/incorrect BibTeX metadata you notice.

- `mucherino2012discretizable` is missing DOI `10.1007/s11590-011-0358-3`.
- `cassioli2015discretization` is missing DOI `10.1016/j.dam.2014.08.035`.
- `abud2018k` is missing DOI `10.1007/s11590-018-1294-2`. Its current `year = 2018`, `volume = 14`, `number = 2`, `pages = 1--14` looks internally inconsistent with the final journal issue; Abud 2024 cites this paper as `Optim. Lett. 14 (2020) 469-482`. Verify final publisher metadata before editing the BibTeX.
- `abud2018k` title should protect `K` in BibTeX if the bibliography style lowercases titles.
- `abud2024impossible` appears complete enough for the fields reviewed here.
