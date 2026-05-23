---
term: Algebraic Exceptional Set
tags: [algebra, geometry]
sources: []
---

## Intuition

An algebraic exceptional set is the nongeneric part of the coordinate parameter
space where a statement that is generically true can fail.

In the DDGP rank-count theory, exceptional sets account for accidental
coincidences such as a vertex lying on a mirror hyperplane even though it is not
part of the mirror clique. Such coincidences can create extra local solutions
that are not predicted by the purely combinatorial rank formula.

## Formal definition

Let the coordinates or distance parameters of a DDGP family be represented by
variables $p$. A property holds outside a proper algebraic exceptional set when
there are nonzero polynomials $f_1,\dots,f_r$ such that the property holds for
all parameters outside

$$
\mathcal E=\{p\mid f_1(p)=\cdots=f_r(p)=0\}.
$$

The set $\mathcal E$ is proper when it is not the whole parameter space. Over
the real numbers, such sets have Lebesgue measure zero whenever the ambient
parameter space is full-dimensional.
