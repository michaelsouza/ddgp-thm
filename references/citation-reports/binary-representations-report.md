# Binary Representation Citation Recommendation

## 1. Papers Reviewed

- Mucherino, "An Analysis on the Degrees of Freedom of Binary Representations for Solutions to Discretizable Distance Geometry Problems" (`mucherino2020analysis`).
- MacNeil and Bodur, "Constraint programming approaches for the discretizable molecular distance geometry problem" (`macneil2022constraint`).
- de Salles Neto, Lavor, and Lodwick, "A note on the Cayley-Menger determinant and the Molecular Distance Geometry Problem" (`de2021note`).
- Baez Sanchez and Lavor, "On the estimation of unknown distances for a class of Euclidean distance matrix completion problems with interval data" (`sanchez2020estimation`).

## 2. Citation-Supported Claims We Can Safely Make

- Exact-distance DDGP/BP search domains can be represented as binary trees, and a solution path can be encoded by binary branch variables. Suggested key: `mucherino2020analysis`; `de2021note` also supports the DMDGP-specific binary-tree statement.
- In DMDGP/CTOP, the key extra structure is consecutivity: the first `K` vertices form a clique and each later vertex is adjacent to the `K` immediately preceding vertices, equivalently yielding overlapping `(K+1)`-cliques. Suggested key: `macneil2022constraint`.
- For any exact-distance DDGP binary path representation, the initial clique bits and the layer `K+1` symmetry can be ignored, leaving the standard `|V|-K-1` branch degrees of freedom. Suggested key: `mucherino2020analysis`.
- Under consecutivity, pruning distances can delimit consecutive subinstances whose dependent binary variables may be replaced by a smaller code block, reducing binary degrees of freedom. Suggested key: `mucherino2020analysis`.
- The current article's graph-derived local branch shifts are adjacent to, but not directly contained in, Mucherino's consecutive-block reduction: Mucherino's note explicitly leaves extension beyond the consecutivity assumption as future work. Suggested key for positioning only: `mucherino2020analysis`.
- Interval/noisy-distance work should be cited only for a scope caveat or future-work sentence, not for the rank-count theorem: `de2021note` supports Cayley-Menger/dual-BP handling of interval data, and `sanchez2020estimation` supports interval EDM-completion estimates and the claim that pruning-phase uncertainty remains challenging.

## 3. Existing Article Citations That Should Remain/Change/Remove

- Keep the existing binary-representation citation to `mucherino2020analysis` in the introduction; it is the strongest reviewed source for branch-bit encodings and reduced degrees of freedom.
- Keep `mucherino2012exploiting` only as a DMDGP symmetry citation, not as the main support for binary degrees of freedom.
- No need to add `macneil2022constraint` to the current binary-code paragraph unless the article adds a sentence about CTOP/DMDGP consecutivity or ordering formulations.
- Do not add `de2021note` or `sanchez2020estimation` to the main rank-count narrative. They are useful if the article adds a short limitation/future-work note saying interval or noisy distances are outside the exact generic setting.
- Avoid wording that says prior binary-representation work already handles non-consecutive DDGP branch-shift rank counts; the reviewed Mucherino paper does not support that.

## 4. Missing/Incorrect BibTeX Metadata

- `mucherino2020analysis`: local notes say the cited version is in *Recent Advances in Computational Optimization*, Springer, pp. 251--255, 2021, with HAL id `hal-03688786`. The current BibTeX says `booktitle={The Workshop on Computational Optimization}` and `year={2020}`; verify and likely update to the published 2021 book/chapter metadata.
- `macneil2022constraint`: DOI missing: `10.1002/net.22068`.
- `de2021note`: DOI missing: `10.1016/j.ins.2020.12.072`. The key `de2021note` also does not match the project's surname-year-title key style; consider `sallesNeto2021note` if keys are normalized later.
- `sanchez2020estimation`: DOI missing: `10.1016/j.laa.2020.01.036`. The first author's surname appears in local notes/references as Baez Sanchez; consider normalizing the author field/key if enforcing project key style, e.g. `baezSanchez2020estimation`.

## 5. Concise Recommendation

Use `mucherino2020analysis` as the central adjacent reference. Add `macneil2022constraint` only for DMDGP/CTOP consecutivity context. Leave `de2021note` and `sanchez2020estimation` out of the main article unless adding a brief exact-vs-interval scope limitation.
