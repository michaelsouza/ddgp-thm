---
term: Canonical Obstruction Map
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The canonical obstruction map assigns a presentation-independent obstruction to
a branch difference.

A branch difference may have several graph-generator presentations. Their raw
labelled violation vectors can differ, so the obstruction is only meaningful
after quotienting out pure presentation labels. The canonical obstruction map is
the resulting well-defined map from branch masks to labelled obstruction
classes.

## Formal definition

Let $M:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{B_P}$ be the mask map, let
$V_F:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{\mathcal L_F}$ be the labelled
violation map, and let

$$
R_F=V_F(\ker M).
$$

The [labelled obstruction quotient](labelled-obstruction-quotient.md) is

$$
Q_F=\mathbb F_2^{\mathcal L_F}/R_F.
$$

For a branch difference $h\in\mathbb F_2^{B_P}$, choose any generator
coefficient vector $\alpha$ such that $M\alpha=h$. The canonical obstruction map
is

$$
\omega_F(h)=V_F\alpha+R_F\in Q_F.
$$

This is independent of the chosen presentation because two presentations of the
same mask differ by an element of $\ker M$.
