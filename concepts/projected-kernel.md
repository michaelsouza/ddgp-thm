---
term: Feasible Branch Shift Space
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The feasible branch shift space is the space of branch-code shifts that admit a
zero-labeled generator presentation.

Applying one of its masks to a feasible local solution preserves every active
edge generically.

## Formal definition

Let $M$ be the mask map and let $V_F$ be the violation map. The feasible branch
shift space is

$$
\mathcal K_F=M(\ker V_F).
$$

Equivalently, using the [canonical obstruction map](canonical-obstruction-map.md),

$$
\mathcal K_F=\ker\omega_F.
$$

When the [labeled presentation property](labeled-presentation-property.md)
holds and $\Xi_F$ is nonempty,

$$
\Xi_F=s^\ast\oplus\mathcal K_F
$$

for any reference solution $s^\ast\in\Xi_F$.

_Avoid_: projected kernel as the article term. It is useful linear-algebra
shorthand, but **feasible branch shift space** names the role of
$\mathcal K_F$ in the solution set.
