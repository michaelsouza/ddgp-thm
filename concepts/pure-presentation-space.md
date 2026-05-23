---
term: Pure Presentation Space
tags: [linear-algebra, geometry, group-theory]
sources: []
---

## Intuition

The pure presentation space records labelled violations that arise only from
changing the generator presentation of the same branch mask.

Such labels are not intrinsic obstructions of a branch difference. They must be
quotiented out to obtain the presentation-independent labelled obstruction.

## Formal definition

Let

$$
M:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{B_P}
$$

be the mask map, and let

$$
V_F:\mathbb F_2^{\mathcal G}\to\mathbb F_2^{\mathcal L_F}
$$

be the labelled violation map. A pure presentation is a generator combination
$\alpha$ with zero total branch mask:

$$
\alpha\in\ker M.
$$

The pure presentation space is

$$
R_F=V_F(\ker M)\subseteq\mathbb F_2^{\mathcal L_F}.
$$

It is the denominator in the [labelled obstruction
quotient](labelled-obstruction-quotient.md)

$$
Q_F=\mathbb F_2^{\mathcal L_F}/R_F.
$$
