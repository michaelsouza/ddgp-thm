# A Simplified Rank-Count Theory for Local Solutions in the Discretizable Distance Geometry Problem

**Author:** Distance Geometry Research Group  
**Status:** Condensed Formal Report  
**Date:** May 2026  

---

## Abstract

This report presents a simplified rank-count theory for the Discretizable Distance Geometry Problem (DDGP). The goal is to count the generic number of local realizations satisfying a prescribed set of active distance constraints without enumerating the Branch-and-Prune tree. The method replaces branch enumeration with linear algebra over $\mathbb F_2$.

For a fixed active edge set $F$, we construct graph-based reflection generators, collect their branch masks in a matrix $M$, and record their labeled active-edge incompatibilities in a matrix $V$. Under genericity assumptions, and assuming that $F$ contains the local base edges, the number of valid local branch codes is

$$|\Xi_F|
=
2^{\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}-\operatorname{rank}(V)}.$$

All ranks are computed over $\mathbb F_2$.

---

## 1. Purpose and Context

The DDGP is a structured form of the Distance Geometry Problem in $\mathbb R^K$. Its vertices are ordered as

$$v_1,\ldots,v_n,$$

and each non-seed vertex $v_i$ is determined by distances to a predecessor set $U_i$ of size $K$. The first $K+1$ vertices form a fixed initial simplex. After these seed vertices, each new vertex normally has two possible positions. Thus, a realization can be represented by a binary branch code.

In the special DMDGP case, predecessor sets are consecutive. This produces interval-reflection symmetries that can be exploited by symmetry-based algorithms. In the general DDGP, predecessor sets need not be consecutive. Reflections therefore no longer act on simple suffixes of the vertex order; they propagate through a dependency graph.

The aim of this theory is to count local valid branch codes for a chosen set of pruning edges $P$, while avoiding exponential enumeration. The central idea is to replace interval symmetries by graph-derived generators and to test their compatibility with active edge constraints by rank computations over $\mathbb F_2$.

---

## 2. Local DDGP Subproblem

Let

$$V_0=\{v_1,\ldots,v_{K+1}\}$$

be the fixed initial simplex. For every non-seed vertex $i>K+1$, the predecessor set $U_i$ has size $K$. The predecessor relation defines a directed acyclic graph: there is an arc $u\to i$ whenever $u\in U_i$.

The **fixing set** of a vertex $i$ is the set of vertices needed to determine it:

$$\operatorname{Fix}_U(i)
=
\{i\}\cup \bigcup_{u\in U_i}\operatorname{Fix}_U(u).$$

The **dependency cone** of a vertex $i$ is the set of vertices whose coordinates depend on the branch choice at $i$:

$$\operatorname{Cone}_U(i)
=
\{i\}\cup \bigcup_{w:\,i\in U_w}\operatorname{Cone}_U(w).$$

For a set of pruning edges $P$, define the local vertex set

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
\{s\in \mathbb F_2^{B_P}: s \text{ realizes every edge in }F\}.$$

The task is to compute $|\Xi_F|$.

---

## 3. Graph-Based Reflection Generators

A generator represents a branch transformation together with the mirror clique that explains its geometric action. Each generator is a pair

$$g=(m_g,C_g),$$

where $m_g\in\mathbb F_2^{B_P}$ is a branch mask and $C_g$ is a predecessor clique. The mask tells which local branch bits are toggled; the clique tells which affine mirror is used.

Let $\mathcal G$ be the set of generators, and write

$$A=\mathbb F_2^{\mathcal G}$$

for the generator coefficient space. An element $\alpha\in A$ is a binary combination of generators.

### 3.1 Cone Generators

For each decision vertex $q\in B_P$, the cone generator is

$$g_q=(m_q,U_q),$$

with support

$$\operatorname{supp}(m_q)=B_P\cap \operatorname{Cone}_U(q).$$

This represents changing the branch choice at $q$ and propagating its effect to all later vertices depending on $q$. Ordered by the DDGP order, the cone masks form a triangular system with ones on the diagonal. Therefore, they span the full local branch space $\mathbb F_2^{B_P}$.

### 3.2 Base Generators

For a predecessor clique $C$, define

$$\operatorname{Gen}_P(C)
=
\{w\in B_P: U_w=C\}.$$

This is the set of local decision vertices generated from the same predecessor clique $C$.

The base generator associated with $C$ is

$$g_C=(m_C,C),$$

where the support of $m_C$ is the union of the cone supports rooted at vertices generated from $C$:

$$\operatorname{supp}(m_C)
=
\bigcup_{w\in \operatorname{Gen}_P(C)}
\operatorname{supp}(m_w).$$

Base generators are useful because two transformations with the same branch mask may have different geometric meaning if they use different mirror cliques.

### 3.3 Generator Mask Matrix

Let

$$\mathcal G=\{g^1,\ldots,g^m\}.$$

The generator mask matrix is the linear map

$$M:A\to \mathbb F_2^{B_P},$$

whose $j$-th column is $m_{g^j}$. Thus

$$M\alpha
=
\bigoplus_{g^j:\,\alpha_j=1}m_{g^j}.$$

The vector $M\alpha$ is the branch mask produced by the generator combination $\alpha$.

---

## 4. Labeled Violations

Fix the active edge set $F$. Let $g=(m_g,C_g)$, and write

$$S_g=\operatorname{supp}(m_g).$$

An active edge $e=\{a,b\}\in F$ crosses the support boundary of $g$ if exactly one endpoint lies in $S_g$:

$$|e\cap S_g|=1.$$

If this happens, one endpoint is moved by the reflection and the other is fixed. The edge length is preserved only when the fixed endpoint lies on the reflection mirror, which is ensured when it belongs to $C_g$. Therefore, $g$ violates $e$ if

$$|e\cap S_g|=1
\quad\text{and}\quad
 e\setminus S_g\not\subseteq C_g.$$

Since $|e\cap S_g|=1$, the set $e\setminus S_g$ is a singleton. Thus the second condition simply says that the fixed endpoint is not in $C_g$.

Violations are labeled by both the edge and the mirror clique. Let

$$\mathcal L
=
\{(e,C_g): g\in\mathcal G,\ e\in F,\ g\text{ violates }e\}.$$

The labeled violation matrix is the linear map

$$V:A\to \mathbb F_2^{\mathcal L}.$$

Rows are indexed by labeled pairs $\ell_i=(e_i,C_i)\in\mathcal L$, and columns are indexed by generators $g^j$. Its entries are

$$V_{ij}=1
\quad\Longleftrightarrow\quad
C_i=C_{g^j}\text{ and }g^j\text{ violates }e_i.$$

Equivalently,

$$V_{ij}=1
\quad\Longleftrightarrow\quad
C_i=C_{g^j},\quad
|e_i\cap S_{g^j}|=1,\quad
 e_i\setminus S_{g^j}\not\subseteq C_{g^j}.$$

Thus $V\alpha$ is the net labeled violation pattern of the generator combination $\alpha$. The same active edge may be crossed by transformations using different mirrors; such violations are kept in different rows and therefore cannot cancel unless their mirror labels agree.

---

## 5. Feasible Branch Shifts and Presentation Redundancy

The branch transformations that preserve all active edges are represented by generator combinations with zero net labeled violation. The feasible branch shift space is

$$\mathcal K_F
=
M(\ker V).$$

Thus, $\mathcal K_F$ consists of branch shifts that can be produced by zero-violation generator combinations.

It is important not to count $\ker V$ itself. Elements of $\ker V$ are generator presentations, not branch masks. Different generator combinations may produce the same branch mask. If

$$M\alpha=M\beta,$$

then

$$\alpha+\beta\in\ker M.$$

So the redundancy is caused by generator combinations that have no net effect on the branch code. The object to count is therefore the image

$$M(\ker V),$$

not the set of presentations $\ker V$.

For the intrinsic obstruction theory, one may quotient out pure presentation artefacts by

$$Q_F
=
\mathbb F_2^{\mathcal L}/V(\ker M).$$

This gives a canonical obstruction map

$$\omega_F:\mathbb F_2^{B_P}\to Q_F,$$

whose kernel is exactly

$$\ker\omega_F=\mathcal K_F.$$

For the rank count, however, it is enough to compute the dimension of $M(\ker V)$ directly.

---

## 6. Main Results

### Branch-Shift Rank Lemma

For the feasible branch shift space

$$\mathcal K_F=M(\ker V),$$

the dimension is

$$\dim\mathcal K_F
=
\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}
-
\operatorname{rank}(V).$$

Proof. Let

$$A=\mathbb F_2^{\mathcal G},
\qquad
m=\dim A=|\mathcal G|.$$

Since $\mathcal K_F=M(\ker V)$, restrict $M$ to $\ker V$:

$$M|_{\ker V}:\ker V\to \mathbb F_2^{B_P}.$$

By rank-nullity applied to this restricted map,

$$\dim\mathcal K_F
=
\dim\ker V
-
\dim(\ker V\cap\ker M).$$

Now

$$\dim\ker V=m-\operatorname{rank}(V),$$

and

$$\ker V\cap\ker M
=
\ker\begin{bmatrix}M\\ V\end{bmatrix}.$$

Therefore,

$$\dim(\ker V\cap\ker M)
=
m-
\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}.$$

Substituting gives

$$\begin{aligned}
\dim\mathcal K_F
&=
\left(m-\operatorname{rank}(V)\right)
-
\left(m-
\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}
\right) \\
&=
\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}
-
\operatorname{rank}(V).
\end{aligned}$$

This proves the rank formula for $\dim\mathcal K_F$.

### Feasible Branch Shift Theorem

Assume:

1. $K\ge 2$;
2. $E_0[L_P]\subseteq F$;
3. $\Xi_F\neq\varnothing$;
4. the DDGP instance is generic, so accidental algebraic degeneracies are excluded.

Then, for any reference solution $s^\ast\in\Xi_F$,

$$\Xi_F
=
s^\ast\oplus \mathcal K_F
=
s^\ast\oplus M(\ker V).$$

Thus the feasible local branch codes form an affine coset of the feasible branch shift space.

### Rank-Count Corollary

Under the assumptions of the Feasible Branch Shift Theorem,

$$|\Xi_F|=|\mathcal K_F|=2^{\dim \mathcal K_F}.$$

By the Branch-Shift Rank Lemma, the local solution count is

$$|\Xi_F|
=
2^{
\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}
-
\operatorname{rank}(V)
}.$$

---

## 7. Counting Algorithm

The rank-count algorithm is entirely combinatorial and works over $\mathbb F_2$.

```text
Input:
  DDGP order, predecessor sets U_i, base edges E_0, pruning edges P

1. Compute Fix_U(i) for all vertices.
2. Build L_P from the fixing sets of the endpoints of P.
3. Set B_P = L_P \ V_0.
4. Set F = E_0[L_P] union P.
5. Build all cone generators g_q for q in B_P.
6. Build all base generators g_C for distinct predecessor cliques C.
7. Construct M from the generator masks.
8. Construct V by applying the labeled violation rule to all active edges.
9. Compute ranks over GF(2):
       r = rank([M; V]) - rank(V).
10. Return 2^r.
```

The method avoids enumerating the $2^{|B_P|}$ local branch strings.

---

## 8. Example Summary

In the running $K=2$ example from the full report, the local decision set is

$$B_P=\{v_4,v_5,v_6\}.$$

The generator mask matrix is

$$M=
\begin{bmatrix}
1&0&0&1&0\\
0&1&0&1&0\\
1&0&1&1&1
\end{bmatrix},$$

and the labeled violation matrix is

$$V=
\begin{bmatrix}
1&1&0&0&0\\
0&0&1&0&1
\end{bmatrix}.$$

The ranks are

$$\operatorname{rank}
\begin{bmatrix}
M\\V
\end{bmatrix}
=3,
\qquad
\operatorname{rank}(V)=2.$$

Thus,

$$|\Xi_F|=2^{3-2}=2.$$

The two feasible local solutions differ by a simultaneous reflection of the relevant local vertices through the appropriate predecessor mirror.

---

## 9. Validation and Scope

The example above is intended as a small reproducible check of the construction: the generator masks, labeled violation matrix, ranks, and final count can all be inspected directly.

The theory is generic. In non-generic configurations, the actual number of solutions can be larger than the generic rank count. This happens when accidental geometric coincidences cause an apparently violating edge to be preserved, for example because a point lies unexpectedly on a reflection mirror. Such cases form a proper algebraic exceptional set.

The case $K=1$ requires separate treatment because the bracket-separation argument used for $K\ge 2$ degenerates on the line.

---

## 10. Conclusion

The rank-count theory provides a concise algebraic method for counting generic local DDGP solutions. It replaces exponential branch enumeration with a finite $\mathbb F_2$ rank computation. The main construction has three ingredients:

1. graph-based branch-mask generators;
2. labeled active-edge violations;
3. the feasible branch shift space $\mathcal K_F=M(\ker V)$.

Under the stated genericity assumptions, the feasible local branch codes form a coset of this shift space, and the number of local solutions is

$$|\Xi_F|
=
2^{
\operatorname{rank}\begin{bmatrix}M\\ V\end{bmatrix}
-
\operatorname{rank}(V)
}.$$

This gives a practical and theoretically justified counting method for general DDGP instances with non-consecutive predecessor structure.

---

## References

1. Lavor, C., Liberti, L., Maculan, N., & Mucherino, A. (2012). *The discretizable distance geometry problem*. Computational Optimization and Applications, 52, 115–146.
2. Liberti, L., Lavor, C., Maculan, N., & Mucherino, A. (2011). *Double discretization in distance geometry*. arXiv:1111.2225.
3. Gonçalves, D. S., & Mucherino, A. (2021). *A new symmetry-based build-up algorithm for the discretizable molecular distance geometry problem*. Journal of Global Optimization, 81, 1011–1033.
