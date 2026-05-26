---
term: Labeled Violation Kernel
tags: [geometry, group-theory, combinatorial-optimization]
sources: []
---

## Intuition

The labeled violation kernel is a linear-algebraic object used to combine local graph-mask symmetries in the presence of active edge constraints.

For one edge, it is often enough to keep graph masks that are individually [mirror-compatible](mirror-compatible-graph-mask.md) with that edge. With multiple active edges, this is too restrictive: two masks may each violate the same edge, but their composition may cancel the violation.

The cancellation is not purely by edge parity. Experiments show that the violation should be labeled by both:

- the active edge being violated;
- the mirror clique of the graph mask causing the violation.

Thus, two violations cancel only when they have the same labeled type.

## Formal definition

Let $\mathcal{G}$ be a family of graph-derived masks

$$
g=(m_g,C_g),
$$

where $m_g\in\mathbb{F}_2^{B}$ is the binary mask and $C_g$ is its mirror clique.

Let $F$ be a set of active edges. Define the labeled violation set

$$
\mathcal{L}
=
\{(e,C_g) \mid g\in\mathcal{G},\ e\in F,\ g \text{ is not mirror-compatible with } e\}.
$$

The violation signature of a generator $g$ is the vector

$$
\nu(g)\in\mathbb{F}_2^{\mathcal{L}},
$$

where

$$
\nu(g)_{(e,C)}
=
1
\quad\Longleftrightarrow\quad
C=C_g
\text{ and }
g \text{ is not mirror-compatible with } e.
$$

Extend $\nu$ linearly to combinations of generators. The labeled violation kernel is

$$
\ker \nu
=
\{\alpha\in\mathbb{F}_2^{\mathcal{G}} \mid \nu(\alpha)=0\}.
$$

The corresponding mask space is the projection

$$
\mathcal{K}_F
=
\left\{
\bigoplus_{g\in\mathcal{G}}\alpha_g m_g
\ \middle|\
\alpha\in\ker\nu
\right\}
\subseteq \mathbb{F}_2^B.
$$

## Rank identity

Let $M$ be the matrix whose columns are the support masks $m_g$, and let $V$ be the matrix whose columns are the labeled violation vectors $\nu(g)$.

The rank of the feasible branch shift space can be computed without enumerating
the kernel:

$$
\operatorname{rank}(\mathcal{K}_F)
=
\operatorname{rank}
\begin{bmatrix}
M\\
V
\end{bmatrix}
-
\operatorname{rank}(V).
$$

This follows by viewing the column span of $\begin{bmatrix}M\\V\end{bmatrix}$ and taking its intersection with the subspace whose violation coordinates are zero.

This rank identity is the algebraic core of the [DDGP labeled-violation rank count](ddgp-labeled-violation-rank-count.md). The geometric part requires the [labeled presentation property](labeled-presentation-property.md): violations with different labels should not cancel outside exceptional coordinate sets. This is the same genericity issue described by [generic mirror separation](generic-mirror-separation.md).

When the same branch mask has multiple generator presentations, the invariant
obstruction lives in the [labeled obstruction quotient](labeled-obstruction-quotient.md), not in the raw labeled violation space.

## Counting role

For generic DDGP instances with an active local edge set $F$, the labeled
presentation theorem gives

$$
|\Xi_F|
=
2^{\operatorname{rank}(\mathcal{K}_F)}.
$$

For local pruning-count experiments, $F=E_0[L_P]\cup P$. This extends the one-pruning-edge mirror-compatible mask rule to base and pruning edges simultaneously. The corresponding algorithmic form is the [DDGP labeled-violation rank count](ddgp-labeled-violation-rank-count.md).
