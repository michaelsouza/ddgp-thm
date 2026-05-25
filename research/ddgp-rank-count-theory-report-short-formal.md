# A Simplified Rank-Count Theory for Local Solutions in the Discretizable Distance Geometry Problem

**Author:** Distance Geometry Research Group  
**Status:** Condensed Formal Report  
**Date:** May 2026  

---

## Abstract

This report presents a simplified formulation of a rank-count theory for the Discretizable Distance Geometry Problem (DDGP). The objective is to count the generic number of local realizations satisfying a set of active distance constraints without enumerating the full Branch-and-Prune tree. The method replaces branch enumeration with linear algebra over $\mathbb{F}_2$. It constructs graph-based reflection generators, records their active-edge incompatibilities in a labelled violation matrix, and obtains the local solution count from the rank formula

$$|\Xi_F|
=
2^{\operatorname{rank}\begin{bmatrix}M\\ V_F\end{bmatrix}
-
\operatorname{rank}(V_F)}.$$

Here $M$ records branch-mask generators, and $V_F$ records labelled active-edge violations. Under genericity assumptions and assuming the active edge set contains the local base edges, this formula gives the number of valid local branch codes.

---

## 1. Purpose and Context

The DDGP is a structured form of the Distance Geometry Problem in $\mathbb{R}^K$. Its vertices are ordered as

$$v_1,\ldots,v_n,$$

and each non-seed vertex $v_i$ is determined by distances to a predecessor set $U_i$ of size $K$. The first $K+1$ vertices form a fixed initial simplex. After these seed vertices, each new vertex normally has two possible positions. Thus, a realization can be represented by a binary branch code.

In the special DMDGP case, predecessor sets are consecutive. This creates a simple interval-reflection structure that can be exploited by symmetry-based algorithms. In the general DDGP, predecessor sets need not be consecutive. Consequently, reflections no longer act on simple suffixes of the vertex order. They propagate through a dependency graph.

The aim of this theory is to count local valid branch codes for a chosen set of pruning edges $P$, while avoiding exponential enumeration. The central idea is to replace interval symmetries by graph-derived generators and to test their compatibility with the active edge constraints by rank computations over $\mathbb{F}_2$.

---

## 2. Local DDGP Subproblem

Let $V_0=\{v_1,\ldots,v_{K+1}\}$ be the fixed initial simplex. For every non-seed vertex $i>K+1$, the predecessor set $U_i$ has size $K$. The predecessor relation defines a directed acyclic graph: there is an arc $u\to i$ whenever $u\in U_i$.

Two graph sets are used throughout the construction.

First, the **fixing set** of a vertex $i$ is the set of vertices needed to determine it:

$$\operatorname{Fix}_U(i)
=
\{i\}\cup \bigcup_{u\in U_i}\operatorname{Fix}_U(u).$$

Second, the **dependency cone** of a vertex $i$ is the set of vertices whose coordinates depend on the branch choice at $i$:

$$\operatorname{Cone}_U(i)
=
\{i\}\cup \bigcup_{w:\,i\in U_w}\operatorname{Cone}_U(w).$$

For a set of pruning edges $P$, the local vertex set is

$$L_P
=
\bigcup_{\{u,v\}\in P}
\left(\operatorname{Fix}_U(u)\cup \operatorname{Fix}_U(v)\right).$$

The local binary decision set is

$$B_P=L_P\setminus V_0.$$

The active edge set is denoted by $F$. In the main theorem, $F$ must contain the local base edges $E_0[L_P]$. In applications, one usually takes

$$F=E_0[L_P]\cup P.$$

The local solution code is

$$\Xi_F
=
\{s\in \mathbb{F}_2^{B_P}:\ s \text{ realizes every edge in }F\}.$$

The task is to compute $|\Xi_F|$.

---

## 3. Graph-Based Reflection Generators

A generator represents a branch transformation together with the mirror clique that explains its geometric action. Each generator is a pair

$$g=(m_g,C_g),$$

where $m_g\in\mathbb{F}_2^{B_P}$ is a branch mask and $C_g$ is a predecessor clique. The support of $m_g$ is the set of branch bits toggled by $g$.

Two types of generators are used.

### 3.1 Dependency-Cone Generators

For each decision vertex $q\in B_P$, the cone generator is

$$g_q=(m_q,U_q),$$

with support

$$\operatorname{supp}(m_q)=B_P\cap \operatorname{Cone}_U(q).$$

This represents changing the branch choice at $q$ and propagating its effect to all later vertices depending on $q$. Ordered by the DDGP order, the cone masks form a triangular system with ones on the diagonal. Therefore, they span the full local branch space $\mathbb{F}_2^{B_P}$.

### 3.2 Base-Clique Closure Generators

For a predecessor clique $C$, let

$$\operatorname{Gen}_P(C)
=
\{w\in L_P\setminus V_0:\ U_w=C\}.$$

The closure generator associated with $C$ has mirror $C$ and support equal to the union of the cone supports rooted at vertices generated from $C$:

$$\operatorname{supp}(m_C)
=
\bigcup_{w\in \operatorname{Gen}_P(C)}
\operatorname{supp}(m_w).$$

Closure generators are useful because two transformations with the same branch mask may have different geometric meaning if they use different mirror cliques.

### 3.3 Generator Mask Matrix

Let $\mathcal{G}$ be the set of all cone and closure generators. The generator mask matrix is

$$M
=
\begin{bmatrix}
| & | & & |\\
m_{g_1} & m_{g_2} & \cdots & m_{g_m}\\
| & | & & |
\end{bmatrix},$$

a binary matrix with one column per generator. It maps generator combinations to branch masks:

$$M\alpha
=
\bigoplus_{g\in \operatorname{supp}(\alpha)}m_g.$$

---

## 4. Labelled Violations

A generator may fail to preserve an active edge. Let $g=(m_g,C_g)$, and let $S_g=\operatorname{supp}(m_g)$. An active edge $e=\{a,b\}\in F$ crosses the support boundary of $g$ if exactly one endpoint lies in $S_g$:

$$|e\cap S_g|=1.$$

If this happens, one endpoint is moved by the reflection and the other is fixed. The edge length is preserved only when the fixed endpoint lies on the reflection mirror, which is ensured when it belongs to $C_g$. Therefore, $g$ violates $e$ if

$$|e\cap S_g|=1
\quad\text{and}\quad
e\setminus S_g\not\subseteq C_g.$$

Violations are labelled by both the edge and the mirror clique:

$$(e,C_g).$$

This distinction is important. The same edge may be crossed by transformations using different mirrors, and such violations should not cancel unless their mirror labels agree.

The labelled violation matrix $V_F$ has one column per generator and one row per labelled violation. Its entry is one precisely when the generator produces that labelled violation. Thus,

$$V_F\alpha$$

is the net labelled violation pattern of the generator combination $\alpha$.

---

## 5. Stabilizer Space and Obstruction

Some generator combinations preserve every active edge because their labelled violations cancel. The stabilizer space is

$$\mathcal{K}_F
=
M(\ker V_F).$$

Thus, $\mathcal{K}_F$ consists of branch masks that can be produced by generator combinations with zero net labelled violation.

There is one technical issue: the same branch mask may have more than one generator presentation. If $\alpha$ and $\beta$ have the same mask, then $\alpha+\beta\in\ker M$. Any violations produced by elements of $\ker M$ are artefacts of presentation rather than intrinsic obstructions. These are removed by the quotient

$$Q_F
=
\mathbb{F}_2^{\mathcal{L}_F}/V_F(\ker M).$$

This quotient defines a canonical obstruction map

$$\omega_F:\mathbb{F}_2^{B_P}\to Q_F.$$

A branch mask has zero obstruction exactly when it belongs to the stabilizer space:

$$\ker \omega_F=\mathcal{K}_F.$$

Informally, $\omega_F(h)=0$ means that the branch change $h$ has no genuine labelled obstruction after presentation artefacts are removed.

---

## 6. Main Rank-Count Theorem

Assume:

1. $K\ge 2$;
2. $E_0[L_P]\subseteq F$;
3. $\Xi_F\neq \varnothing$;
4. the DDGP instance is generic, meaning that accidental algebraic degeneracies are excluded.

Then, for any feasible local code $s_0\in\Xi_F$,

$$\Xi_F
=
s_0\oplus M(\ker V_F).$$

Therefore, the local solution code is a coset of the stabilizer space. Its size is

$$|\Xi_F|
=
2^{\dim \mathcal{K}_F}.$$

The dimension of $\mathcal{K}_F$ is computed by rank-nullity. Define

$$T(\alpha)=(M\alpha,V_F\alpha).$$

Then

$$\dim \mathcal{K}_F
=
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F).$$

Hence,

$$|\Xi_F|
=
2^{
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F)
}.$$

This is the main counting formula.

---

## 7. Simplified Proof Idea

The proof has four main components.

First, each non-seed vertex has two possible positions, and these positions are reflections through the affine span of its predecessor clique. This gives the geometric basis for the generator model.

Second, a generator preserves an active edge unless it moves exactly one endpoint while the fixed endpoint is not on the mirror. This is exactly the labelled violation rule.

Third, if a generator combination has zero net labelled violation, then its transformations can be grouped by mirror clique into admissible reflection blocks. These blocks preserve all base edges and all active edges. Therefore, any mask in $M(\ker V_F)$ maps a feasible local code to another feasible local code.

Fourth, if a branch mask has a nonzero canonical obstruction, then at least one active-edge length changes generically. The proof uses algebraic separation of labelled obstruction terms. For $K\ge 2$, distinct labelled obstructions have distinct primitive bracket signatures, so they cannot cancel identically outside a proper algebraic exceptional set.

Together, these facts show that the feasible branch differences are exactly $M(\ker V_F)$, which gives the coset structure and the rank formula.

---

## 8. Counting Algorithm

The rank-count algorithm is entirely combinatorial and works over $\mathbb{F}_2$.

```text
Input:
  DDGP order, predecessor sets U_i, base edges E_0, pruning edges P

1. Compute Fix_U(i) for all vertices.
2. Build L_P from the fixing sets of the endpoints of P.
3. Set B_P = L_P \ V_0.
4. Set F = E_0[L_P] union P.
5. Build all cone generators g_q for q in B_P.
6. Build all closure generators g_C for distinct predecessor cliques C.
7. Construct M from the generator masks.
8. Construct V_F by applying the labelled violation rule to all active edges.
9. Compute ranks over GF(2):
       r = rank([M; V_F]) - rank(V_F).
10. Return 2^r.
```

The method avoids enumerating the $2^{|B_P|}$ local branch strings.

---

## 9. Example Summary

In the running $K=2$ example from the full report, the local decision set is

$$B_P=\{v_4,v_5,v_6\}.$$

The generator mask matrix is

$$M=
\begin{bmatrix}
1&0&0&1&0\\
0&1&0&1&0\\
1&0&1&1&1
\end{bmatrix},$$

and the labelled violation matrix is

$$V_F=
\begin{bmatrix}
1&1&0&0&0\\
0&0&1&0&1
\end{bmatrix}.$$

The ranks are

$$\operatorname{rank}
\begin{bmatrix}
M\\V_F
\end{bmatrix}
=3,
\qquad
\operatorname{rank}(V_F)=2.$$

Thus,

$$|\Xi_F|=2^{3-2}=2.$$

The two feasible local solutions differ by a simultaneous reflection of the relevant local vertices through the appropriate predecessor mirror.

---

## 10. Validation and Scope

The formula was tested against exact coordinate enumeration in dimensions $K=2$ and $K=3$. Across the stored generic test suite, it matched both exact local base counts and exact valid counts in all reported cases:

$$2225/2225.$$

The theory is generic. In non-generic configurations, the actual number of solutions can be larger than the predicted count. This happens when accidental geometric coincidences cause an apparently violating edge to be preserved, for example because a point lies unexpectedly on a reflection mirror. Such cases form a proper algebraic exceptional set.

The case $K=1$ requires separate treatment because the bracket-separation argument used for $K\ge 2$ degenerates on the line.

---

## 11. Conclusion

The rank-count theory provides a concise algebraic method for counting generic local DDGP solutions. It replaces exponential branch enumeration with a finite $\mathbb{F}_2$ rank computation. The main construction has three ingredients:

1. graph-based branch-mask generators;
2. labelled active-edge violations;
3. a stabilizer space obtained as $M(\ker V_F)$.

Under the stated genericity assumptions, the feasible local branch codes form a coset of this stabilizer space, and the number of local solutions is given by

$$|\Xi_F|
=
2^{
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F)
}.$$

This gives a practical and theoretically justified counting method for general DDGP instances with non-consecutive predecessor structure.

---

## References

1. Lavor, C., Liberti, L., Maculan, N., & Mucherino, A. (2012). *The discretizable distance geometry problem*. Computational Optimization and Applications, 52, 115–146.
2. Liberti, L., Lavor, C., Maculan, N., & Mucherino, A. (2011). *Double discretization in distance geometry*. arXiv:1111.2225.
3. Gonçalves, D. S., & Mucherino, A. (2021). *A new symmetry-based build-up algorithm for the discretizable molecular distance geometry problem*. Journal of Global Optimization, 81, 1011–1033.
