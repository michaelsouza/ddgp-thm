# Applications/background citation review

## 1. Papers reviewed

- `billinge2018recent`: broad survey of assigned/unassigned distance geometry, rigidity, DDGP with intervals, proteins, nanostructures, and vector geometry. Recommendation: leave uncited unless adding a broad applications sentence.
- `costa2014discretization`: protein side-chain discretization orders and iBP with NMR pruning distances. Recommendation: leave uncited for the current article; cite only for protein side-chain ordering context.
- `cassioli2015algorithm`: iBP for exhaustive sampling of protein conformations satisfying distance constraints, with protein-specific pruning devices. Recommendation: leave uncited unless discussing interval BP applications.
- `costa2017calculating`: constraint interval analysis for uncertainty propagation in iDMDGP from NMR interval distances. Recommendation: leave uncited; it is about interval data, not the article's generic exact-distance rank count.
- `lavor2011computation`: artificial hydrogen backbones that satisfy DMDGP assumptions and reconstruct protein backbones. Recommendation: leave uncited unless discussing hydrogen-only NMR/backbone modeling.
- `lavor2009discretizable`: short early DMDGP/BP formulation for MDGP instances with binary choices. Recommendation: leave uncited because later, fuller sources already cover the article's needed theory.
- `sallaume2013discrete`: BP extension to whole protein structures including side chains. Recommendation: leave uncited unless the article adds whole-protein application context.

## 2. Citation-supported claims we can safely make

- Distance geometry and DDGP methods have molecular/protein and nanostructure applications, including NMR-derived protein distance constraints: `billinge2018recent`; for protein-specific iBP context, `cassioli2015algorithm`.
- NMR protein data are naturally modeled with sparse and sometimes interval-valued distance constraints, motivating iDMDGP/iBP variants: `costa2014discretization`, `cassioli2015algorithm`, `costa2017calculating`.
- Protein side-chain or hydrogen-backbone orderings can be constructed so that molecular instances satisfy DMDGP/iDMDGP discretization assumptions: `lavor2011computation`, `costa2014discretization`, `sallaume2013discrete`.
- The BP/iBP search tree can enumerate feasible conformations in these protein-focused subclasses: `lavor2009discretizable`, `cassioli2015algorithm`, `sallaume2013discrete`.

Action: do not add these claims to the short article unless the introduction intentionally expands from theory motivation into applications/background. The current article is tighter without them.

## 3. Existing article citations

- Keep `cassioli2015discretization` for the current claim about discretization-order theory. Do not replace it with `costa2014discretization`, which is narrower and protein-side-chain-specific.
- Keep the current DMDGP/DDGP theory citations in the introduction and genericity discussion; none of the reviewed application papers is a better primary source for the rank-count theorem's setup.
- No current article citation should be removed because of this review.
- Do not cite `billinge2018recent` merely to show that distance geometry has applications; that would add broad survey context without supporting a technical claim in the article.

## 4. BibTeX metadata notes

- `billinge2018recent`: DOI available in notes but missing from BibTeX: `10.1007/s10479-018-2989-6`.
- `costa2014discretization`: DOI available in notes but missing from BibTeX: `10.1007/s10898-013-0135-1`.
- `cassioli2015algorithm`: DOI available in notes but missing from BibTeX: `10.1186/s12859-015-0451-1`; journal capitalization should be `BMC Bioinformatics`.
- `costa2017calculating`: DOI available in notes but missing from BibTeX: `10.1016/j.ins.2017.06.015`; volume appears to be `415--416`, not just `415`.
- `lavor2011computation`: DOI available in notes but missing from BibTeX: `10.1007/s10898-010-9584-y`.
- `lavor2009discretizable` and `sallaume2013discrete`: no obvious DOI in the provided notes.
