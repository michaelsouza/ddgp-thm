---
term: Projected Kernel
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The projected kernel is the space of branch masks that admit a zero-labelled
generator presentation.

It is the algebraic stabilizer predicted by the rank-count theory: applying one
of its masks to a feasible local solution should preserve every active edge
generically.

## Formal definition

Let $M$ be the mask map and let $V_F$ be the violation map. The projected kernel
is

$$
\mathcal K_F=M(\ker V_F).
$$

Equivalently, using the [canonical obstruction map](canonical-obstruction-map.md),

$$
\mathcal K_F=\ker\omega_F.
$$

When the [labelled presentation property](labelled-presentation-property.md)
holds and $\Xi_F$ is nonempty,

$$
\Xi_F=s_0\oplus\mathcal K_F
$$

for any feasible local code $s_0\in\Xi_F$.
