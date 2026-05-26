# Labeled-Violation Rank Count: Generic Proof Formulation

## Goal

Let \(F\) be a local active edge set in a DDGP instance. For pruning
experiments,

$$
F=E_0[L_P]\cup P,
$$

where \(E_0[L_P]\) is the induced local base edge set and \(P\) is the
pruning-edge set.

The proposed count is

$$
|\Xi_F|
=
2^{r_F},
\qquad
r_F
=
\operatorname{rank}_{\mathbb F_2}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}_{\mathbb F_2}(V_F).
$$

Here:

- \(L_P\) is the local vertex set induced by the fixing sets of the endpoints
  in \(P\);
- \(B_P=L_P\setminus V_0\) is the local binary decision set after fixing the
  first \(K+1\) vertices;
- \(\mathcal G=\{g=(m_g,C_g)\}\) is the graph-generator family, consisting of
  dependency-cone generators and base-clique closure generators;
- \(M\) is the matrix whose columns are the masks \(m_g\);
- \(V_F\) is the labeled violation matrix for the active edge set \(F\);
- labels have the form \((e,C_g)\), where \(e\in F\) and \(C_g\) is the mirror
  clique of \(g\).

This note gives the proof in three layers:

1. an unconditional linear-algebra theorem for \(M(\ker V_F)\);
2. a finite algebraic-exceptional-set theorem showing how algebraic
   edge-preserving differences from a base code become the generic local
   solution code;
3. the geometric labeled-presentation lemma that must hold for the DDGP
   graph generators.

The first two layers are complete. The third layer is now isolated as one
precise mathematical lemma, rather than being hidden inside the experiments.

## Local Setup

Let the DDGP instance be feasible in \(\mathbb R^K\), with order
\(v_1,\dots,v_n\), fixed initial simplex

$$
V_0=\{v_1,\dots,v_{K+1}\},
$$

and predecessor sets \(U_i\) for \(i>K+1\).

For a pruning-edge set \(P\), define

$$
L_P
=
\bigcup_{\{u,v\}\in P}
\left(\operatorname{Fix}_U(u)\cup\operatorname{Fix}_U(v)\right),
\qquad
B_P=L_P\setminus V_0.
$$

For an active edge set \(F\subseteq E[L_P]\), the local solution code is

$$
\Xi_F
=
\{s\in\mathbb F_2^{B_P}\mid s\text{ realizes every edge in }F\}.
$$

The local valid code for pruning set \(P\) is obtained with

$$
F=E_0[L_P]\cup P.
$$

## Graph Generators

The generator family \(\mathcal G\) contains two types.

For every decision vertex \(q\in B_P\), the dependency-cone generator is

$$
g_q=(m_q,U_q),
\qquad
\operatorname{supp}(m_q)
=
B_P\cap \operatorname{Cone}_U(q).
$$

For every predecessor clique \(C\), the base-clique closure generator is

$$
g_C=(m_C,C),
\qquad
\operatorname{supp}(m_C)
=
B_P
\cap
\bigcup_{w:U_w=C}
\operatorname{Cone}_U(w).
$$

Closure generators are not needed to span the branch space. They are needed
because a block of vertices generated from the same mirror clique can have a
single geometric obstruction label even when the same mask is representable as
a sum of individual cones.

## Labeled Violations

For a generator \(g=(m_g,C_g)\), let

$$
S_g=\operatorname{supp}(m_g).
$$

For an active edge \(e=\{a,b\}\in F\), say that \(g\) violates \(e\) when \(e\)
crosses the support boundary and the fixed endpoint is not in the mirror
clique:

$$
|e\cap S_g|=1
\quad\text{and}\quad
e\setminus S_g\not\subseteq C_g.
$$

The labeled violation set is

$$
\mathcal L_F
=
\{(e,C_g)\mid g\in\mathcal G,\ e\in F,\ g\text{ violates }e\}.
$$

The labeled violation vector is

$$
\nu_F(g)\in\mathbb F_2^{\mathcal L_F},
\qquad
\nu_F(g)_{(e,C)}
=1
\Longleftrightarrow
C=C_g
\text{ and }g\text{ violates }e.
$$

Let \(A=\mathbb F_2^{\mathcal G}\) be the coefficient space of graph
generators. Define linear maps

$$
M:A\to\mathbb F_2^{B_P},
\qquad
V_F:A\to\mathbb F_2^{\mathcal L_F},
$$

by linear extension of \(m_g\) and \(\nu_F(g)\).

The predicted stabilizer is

$$
\mathcal K_F=M(\ker V_F).
$$

Equivalently, if

$$
Z_F
=
\operatorname{span}_{\mathbb F_2}
\{(m_g,\nu_F(g))\mid g\in\mathcal G\}
\subseteq
\mathbb F_2^{B_P}\oplus \mathbb F_2^{\mathcal L_F},
$$

then

$$
\mathcal K_F=\{h\in\mathbb F_2^{B_P}\mid (h,0)\in Z_F\}.
$$

This lifted form is the cleanest formulation: a branch difference is allowed
when it admits a graph-generator presentation with zero total labeled
obstruction.

## Canonical Labeled Obstruction Quotient

The labeled vector \(V_F\alpha\) depends on the chosen generator presentation
\(\alpha\). The branch difference

$$
h=M\alpha
$$

does not. Therefore the invariant obstruction attached to \(h\) is not an
element of \(\mathbb F_2^{\mathcal L_F}\) itself, but a coset modulo labeled
relations that change the presentation without changing the branch mask.

Let

$$
R_F=V_F(\ker M)\subseteq \mathbb F_2^{\mathcal L_F}
$$

be the space of pure presentation labels. Define

$$
Q_F=\mathbb F_2^{\mathcal L_F}/R_F.
$$

Because cone masks span the branch space, for every \(h\in\mathbb F_2^{B_P}\)
there exists \(\alpha\in A\) with \(M\alpha=h\). Define

$$
\omega_F(h)=V_F\alpha+R_F\in Q_F.
$$

This is well-defined: if \(M\alpha=M\beta\), then
\(\alpha+\beta\in\ker M\), so \(V_F\alpha+V_F\beta\in R_F\).

The projected kernel is exactly the kernel of this canonical obstruction map:

$$
\mathcal K_F=\ker\omega_F.
$$

Proof.

If \(h\in\mathcal K_F\), then \(h=M\alpha\) for some
\(\alpha\in\ker V_F\). Hence \(\omega_F(h)=0+R_F\).

Conversely, if \(\omega_F(h)=0\), choose \(M\alpha=h\). Since
\(V_F\alpha\in R_F\), there is \(\delta\in\ker M\) with
\(V_F\delta=V_F\alpha\). Then

$$
M(\alpha+\delta)=h,
\qquad
V_F(\alpha+\delta)=0,
$$

so \(h\in M(\ker V_F)=\mathcal K_F\). \(\square\)

Thus the labeled presentation property can be restated without reference to a
particular generator presentation:

$$
h\text{ preserves all active edges generically}
\quad\Longleftrightarrow\quad
\omega_F(h)=0.
$$

## Elementary Geometry

### Lemma 1: Lateration Choices Are Reflections

Let \(C=U_i\) be the predecessor clique of a non-seed vertex \(v_i\). Once the
vertices in \(C\) are fixed, the two possible positions of \(v_i\) are
reflections of each other through

$$
H_C=\operatorname{aff}(C).
$$

Proof.

The two points satisfy the same \(K\) sphere equations centered at the affinely
independent points in \(C\). Subtracting one sphere equation from another gives
\(K-1\) independent linear equations. Their common solution set is a line
orthogonal to \(H_C\). Intersecting that line with one sphere gives two points
symmetric with respect to \(H_C\). \(\square\)

### Lemma 2: Mirror-Compatible Edge Obstruction

Let \(g=(m,C)\), let \(S=\operatorname{supp}(m)\), and consider the ideal
reflection that reflects all vertices in \(S\) through \(H_C\) while leaving
all vertices outside \(S\) fixed. If an edge \(e=\{a,b\}\) does not violate
\(g\), then this ideal reflection preserves the length of \(e\).

Proof.

If \(a,b\notin S\), nothing changes.

If \(a,b\in S\), both endpoints are reflected by the same Euclidean isometry.

If exactly one endpoint is in \(S\), non-violation means the fixed endpoint
lies in \(C\), hence lies on \(H_C\). Reflection preserves distance to every
point of its mirror hyperplane. \(\square\)

This lemma justifies the mirror-compatible rule for one edge. It does not
claim that every graph generator is itself a valid DDGP symmetry. Cone and
closure generators can cross active base edges. Such crossings are precisely
the labeled obstructions recorded in \(V_F\).

### Lemma 3: A Single Violation Changes Length Generically

Let \(e=\{a,b\}\) cross the support of \(g=(m,C)\), and suppose the fixed
endpoint is not in \(C\). Then the ideal reflection associated with \(g\)
changes the length of \(e\) outside a proper algebraic exceptional set.

Proof.

Choose coordinates so that \(H_C=\{x_K=0\}\). Assume \(a\in S\), \(b\notin S\),
and write

$$
a=(a_\parallel,h_a),
\qquad
b=(b_\parallel,h_b).
$$

The reflected point is \(a'=(a_\parallel,-h_a)\). Hence

$$
\|a'-b\|^2-\|a-b\|^2
=
(-h_a-h_b)^2-(h_a-h_b)^2
=
4h_ah_b.
$$

Strict simplex nondegeneracy gives \(h_a\ne 0\). Generic mirror incidence gives
\(h_b\ne 0\) when \(b\notin C\). Therefore the length changes outside the
algebraic set \(h_ah_b=0\). \(\square\)

### Lemma 4: Cone Masks Span the Branch Space

The dependency-cone masks \(\{m_q\mid q\in B_P\}\) span
\(\mathbb F_2^{B_P}\).

Proof.

Order \(B_P\) by the DDGP order. The cone mask \(m_q\) contains \(q\) and only
vertices \(z\ge q\), because dependency arcs always point forward in the
order. Thus the cone-mask matrix is upper triangular with diagonal entries
equal to \(1\). It has full rank \(|B_P|\). \(\square\)

Consequently every branch-code difference \(s\oplus t\) can be represented as
\(M\alpha\) for some \(\alpha\in A\). The issue is not spanning. The issue is
which spanning presentations have zero geometric obstruction.

## Algebraic Parameter Space

Fix the DDGP order, predecessor sets, local vertex set, and active edge set.
Let \(\Omega_U\) be the semialgebraic set of coordinate realizations for which:

- the initial \(K+1\) vertices form a nondegenerate simplex;
- each predecessor clique \(U_i\) is affinely independent;
- every lateration step has two distinct choices;
- no active-edge endpoint lies accidentally on a mirror hyperplane unless this
  is forced by membership in the mirror clique.

The excluded conditions are polynomial equalities or inequalities in the
coordinates. Thus the bad set is algebraic, and \(\Omega_U\) is Zariski dense
whenever it is nonempty.

For a branch code \(s\in\mathbb F_2^{B_P}\), write \(x_i(s;p)\) for the
coordinate of \(v_i\) obtained from parameter point \(p\in\Omega_U\). The
coordinates are algebraic functions on the finite lateration cover.

For an edge \(e=\{a,b\}\), a branch difference \(h\), and a base code \(s\),
define the edge residual

$$
\Delta_{e,s,h}(p)
=
\|x_a(s\oplus h;p)-x_b(s\oplus h;p)\|^2
-
\|x_a(s;p)-x_b(s;p)\|^2.
$$

After clearing the finitely many lateration square roots, each equation
\(\Delta_{e,s,h}=0\) defines an algebraic subset of the parameter space.

### Lemma 5: Finite Exceptional-Set Lemma

For fixed \(F\), \(B_P\), \(U_i\), branch code \(s\), and branch difference
\(h\), suppose \(h\) does not preserve all active edges from \(s\) identically.
Then the set of parameters for which \(h\) nevertheless preserves all active
edges from \(s\) is contained in a proper algebraic subset.

Proof.

By assumption, at least one active edge \(e\in F\) has
\(\Delta_{e,s,h}\) not equal to the zero algebraic function. The zero set of a
nonzero polynomial, after clearing radicals on the finite lateration cover, is
a proper algebraic subset. There are only finitely many branch codes, branch
differences, and active edges, so the union of all such accidental zero sets is
still a proper algebraic subset. \(\square\)

## Algebraic Stabilizer from a Base Code

For a fixed branch code \(s\), define

$$
\mathcal U_F(s)
=
\{h\in\mathbb F_2^{B_P}
\mid
\Delta_{e,s,h}\equiv 0
\text{ for all }e\in F\}.
$$

This is the set of branch differences that preserve every active edge from the
particular code \(s\) as an algebraic identity. Before the labeled
presentation is proved, \(\mathcal U_F(s)\) should not be assumed to be a vector
space or to be independent of \(s\).

### Theorem 1: Generic Solution Codes Are Algebraic-Stabilizer Translates

Assume \(\Xi_F\neq\varnothing\). Outside a proper algebraic exceptional set,
for every \(s_0\in\Xi_F\),

$$
\Xi_F=s_0\oplus \mathcal U_F(s_0).
$$

Proof.

If \(h\in\mathcal U_F(s_0)\), then \(s_0\oplus h\) realizes every active edge
because all active-edge lengths are preserved identically.

Conversely, let \(h\notin\mathcal U_F(s_0)\). By Lemma 5, after excluding the
finite union of proper algebraic accidental-zero sets, at least one active-edge
residual \(\Delta_{e,s_0,h}\) is nonzero at the current parameter point. Hence
\(s_0\oplus h\notin\Xi_F\).

Therefore the feasible local code is exactly the translate
\(s_0\oplus\mathcal U_F(s_0)\). \(\square\)

Thus the generic counting problem is equivalent to the graph-theoretic problem
of proving

$$
\mathcal U_F(s_0)=\mathcal K_F
\quad\text{for every feasible }s_0.
$$

## Labeled Presentation Lemma

The required geometric statement is the following.

### Lemma 6: Generic Labeled Presentation

For a generic DDGP realization with active edge set \(F\), for every feasible
base code \(s_0\) and every branch difference \(h\),

$$
h\in\mathcal U_F(s_0)
\quad\Longleftrightarrow\quad
\omega_F(h)=0.
$$

Equivalently, for every feasible \(s_0\in\Xi_F\),

$$
\Xi_F=s_0\oplus M(\ker V_F).
$$

This is the DDGP analogue of the DMDGP statement that valid solutions form a
coset of the partial-reflection symmetry group.

### What Lemma 6 Must Prove

Because cone masks already span the branch space, any branch difference
\(h\) has at least one graph-generator presentation \(h=M\alpha\). The
canonical object associated with \(h\) is the obstruction coset
\(\omega_F(h)\in Q_F\), not a particular vector \(V_F\alpha\).
Lemma 6 has two directions.

Soundness:

$$
\omega_F(h)=0
\quad\Longrightarrow\quad
h\in\mathcal U_F(s_0)
\text{ for every feasible }s_0.
$$

Equivalently, if \(h\) has at least one zero-labeled presentation, then \(h\)
is an active-edge symmetry from every feasible base code.

Completeness:

$$
h\in\mathcal U_F(s_0)
\quad\Longrightarrow\quad
\omega_F(h)=0.
$$

That is, every algebraic active-edge symmetry from a feasible base code has a
zero-labeled graph presentation.

The only possible failure of completeness is an accidental cancellation between
different labels after passing to the quotient \(Q_F\). Lemma 3 rules out a
single primitive violation. The remaining case is a combination of violations
whose algebraic length changes cancel for special coordinates even though the
obstruction coset \(\omega_F(h)\) is nonzero.

### Edge-Sign Separation Principle

The natural way to prove Lemma 6 is to prove the following more local claim.

For each active edge \(e\), the lateration construction induces a finite set of
edge signs

$$
\eta_{e,C}\in\{\pm1\},
\qquad
(e,C)\in\mathcal L_F,
$$

such that:

1. applying a generator \(g=(m_g,C_g)\) flips exactly the edge signs whose
   labels occur in \(\nu_F(g)\);
2. the squared length of \(e\) is a generic algebraic function of these signs;
3. a nonzero change in any edge-sign character has a nonzero coefficient
   outside a proper algebraic exceptional set.

Under this principle, active-edge preservation is equivalent to preserving all
edge-sign characters. This is exactly the condition \(V_F\alpha=0\). Thus
Lemma 6 follows.

The single-sign coefficient in item 3 is the quantity \(4h_ah_b\) from
Lemma 3. Higher-order cancellations are excluded by the standard algebraic
independence argument: two distinct labeled sign characters define distinct
polynomials on the lateration cover, and a nontrivial linear combination can
vanish identically only if all its character coefficients vanish. The
exceptional set is the finite union of the coefficient-zero hypersurfaces.

This is now the precise mathematical core of the proof. It replaces the vague
phrase "no accidental mirror cancellations" by the concrete requirement that
the labeled edge-sign characters are generically separated.

## Bracket Obstruction Algebra

The edge-sign separation principle can be made more explicit with bracket
coordinates.

For a \(K\)-clique \(C=\{c_1,\dots,c_K\}\) and a vertex \(z\), define the
oriented bracket

$$
[C,z]
=
\det
\begin{bmatrix}
1&1&\cdots&1\\
x_{c_1}&x_{c_2}&\cdots&x_{c_K}&x_z
\end{bmatrix},
$$

where the determinant is taken in homogeneous coordinates. It vanishes exactly
when \(z\in\operatorname{aff}(C)\). Up to a nonzero normalization depending
only on \(C\), the signed height of \(z\) over \(H_C=\operatorname{aff}(C)\) is
\([C,z]\).

For an active edge \(e=\{a,b\}\) and mirror clique \(C\), the primitive
obstruction monomial is

$$
\mu_{e,C}=[C,a][C,b].
$$

If \(a\in C\) or \(b\in C\), then \(\mu_{e,C}=0\), which is exactly the
mirror-compatible crossing case. If neither endpoint lies in \(C\), then
\(\mu_{e,C}\) is generically nonzero.

### Lemma 7: Primitive Bracket Separation

Assume \(K\ge2\). In the coordinate ring of generic \(K\)-dimensional point
configurations, every nonzero obstruction coset in \(Q_F\) has a detecting
primitive bracket combination. Equivalently, if a labeled obstruction class is
nonzero, then some linear functional on labels that annihilates
\(R_F=V_F(\ker M)\) selects a bracket combination that is not the zero
polynomial.

Proof.

Let \(0\ne q\in Q_F\). Choose a linear functional

$$
\varphi:\mathbb F_2^{\mathcal L_F}\to\mathbb F_2
$$

such that \(\varphi(R_F)=0\) and \(\varphi(q)=1\). It is enough to show that
the bracket combination

$$
B_\varphi
=
\sum_{\varphi(e,C)=1}\mu_{e,C}
$$

is not identically zero, where the sum is interpreted over a characteristic
zero coefficient field.

Specialize the generic points to a stretched moment curve

$$
x_i(t)=\left(t^{w_i},t^{2w_i},\dots,t^{Kw_i}\right),
$$

with strictly increasing integer weights \(w_i=N^i\) and \(N\) larger than the
number of vertices. For \(S=\{i_0,\dots,i_K\}\), the bracket \([S]\) specializes
to a Vandermonde product

$$
\prod_{p<q}\left(t^{w_q}-t^{w_p}\right),
$$

up to sign. Its leading exponent is

$$
E(S)=\sum_{p<q}\max(w_p,w_q).
$$

For a label \((e=\{a,b\},C)\), the primitive obstruction has leading exponent

$$
E(C\cup\{a\})+E(C\cup\{b\}).
$$

Because the weights are powers of \(N\), this exponent encodes the multiset of
unordered vertex pairs appearing in the two brackets. For \(K\ge2\), that
pair-multiset determines the label: the pairs internal to \(C\) are exactly the
pairs appearing with multiplicity two, and the two endpoints are recovered as
the vertices joined to all of \(C\) by multiplicity-one pairs. Thus distinct
labels have distinct leading exponents under this specialization.

Therefore every nonempty sum of primitive obstruction monomials has a unique
largest leading exponent and cannot vanish identically. In particular,
\(B_\varphi\ne0\). \(\square\)

### Lemma 8: Nonzero Labeled Obstruction Gives a Nonzero Edge Residual

Let \(h\in\mathbb F_2^{B_P}\). If \(\omega_F(h)\ne0\), then for generic
coordinates there exists \(e\in F\) such that

$$
\Delta_{e,s_0,h}\not\equiv0.
$$

Proof.

Choose a generator presentation \(h=M\alpha\). Group the selected generators by
their mirror clique. For a fixed active edge \(e\) and mirror clique \(C\), the
parity of non-compatible crossings is the \((e,C)\)-coordinate of
\(V_F\alpha\). Changing the presentation modifies this vector only by an
element of \(R_F\), so the invariant part is \(\omega_F(h)\).

Assume \(\omega_F(h)\ne0\). By Lemma 7, some active edge \(e\) has a nonzero
primitive bracket monomial \(\mu_{e,C}\) in the obstruction expansion. The
coefficient of this monomial in the edge residual is the single-reflection
coefficient from Lemma 3, namely a nonzero scalar multiple of
\(4h_a h_b\). Since \(\mu_{e,C}\) is not cancelled by any other label in the
same obstruction coset, this coefficient cannot vanish identically.

Thus \(\Delta_{e,s_0,h}\) is a nonzero algebraic function. Its zero set is a
proper algebraic subset. \(\square\)

### Lemma 9: Zero Labeled Obstruction Preserves Active Edges

If \(\omega_F(h)=0\), then \(h\in\mathcal U_F(s_0)\) for every feasible
branch code \(s_0\).

Proof.

The condition \(\omega_F(h)=0\) is equivalent to \(h\in M(\ker V_F)\). Choose
\(\alpha\in\ker V_F\) with \(M\alpha=h\). Group the selected generators in
\(\alpha\) by mirror clique \(C\), and let \(m_C^\alpha\) be the XOR of their
support masks.

For every active edge \(e\) and every mirror clique \(C\), the labeled
violation parity is zero. Hence the block mask \((m_C^\alpha,C)\) is
mirror-compatible with \(e\): either both endpoints are toggled by the block,
neither endpoint is toggled, or the only fixed endpoint of a crossing lies in
\(C\). By Lemma 2, reflecting the block through the current hyperplane
\(\operatorname{aff}(C)\) preserves the length of \(e\).

Apply these block reflections one mirror clique at a time, always using the
current position of the mirror clique. Each step preserves every active edge,
so their composition preserves every active edge. The resulting branch-code
difference is the XOR of the block masks, namely

$$
\bigoplus_C m_C^\alpha=M\alpha=h.
$$

Therefore \(s_0\oplus h\) realizes all active edges whenever \(s_0\) does, and
\(h\in\mathcal U_F(s_0)\). \(\square\)

Lemmas 8 and 9 prove Lemma 6 outside the finite algebraic exceptional set where
one of the primitive bracket coefficients vanishes.

## Field-Theoretic Reformulation of Lemma 6

There is an equivalent way to read the obstruction theorem that may be cleaner
for a final paper proof.

Work over a generic coordinate/distance parameter field \(\mathbb K_0\) for the
fixed initial simplex and the prescribed DDGP predecessor graph. Each
lateration step adjoins one quadratic sign choice. Thus, outside the usual
discriminant-zero exceptional set, the full discretization construction gives a
tower or, more canonically, a finite iterated quadratic cover

$$
\mathbb K_0
\subset
\mathbb K_{K+2}
\subset
\cdots
\subset
\mathbb K_n,
$$

where each non-seed vertex contributes one binary lateration sign. The full
splitting field of an iterated quadratic tower need not have an abelian Galois
group. What is canonical for BP is the branch-code torsor
\(\mathbb F_2^{B_P}\), not a priori an abelian group of Euclidean motions.

For a decision vertex \(q\), the nontrivial automorphism of the quadratic step
at \(q\), when transported along the BP cover, reflects \(v_q\) across
\(H_{U_q}\) and propagates to the dependency descendants of \(q\). On branch
codes, this gives the dependency-cone mask. Base-clique closure generators are
the corresponding grouped sign changes for vertices generated by a common
mirror clique.

For an active edge \(e=\{a,b\}\), define the edge-length element

$$
\lambda_e=\|x_a-x_b\|^2\in \mathbb K_n.
$$

A branch difference preserves \(e\) exactly when the two points of the BP cover
labeled by that difference have the same value of \(\lambda_e\). Therefore the
local active-edge stabilizer is the set of branch differences that fix all
functions \(\lambda_e\) on the cover.

The labeled presentation lemma is equivalent to the statement

$$
\operatorname{Stab}_{\mathrm{BP}}(F)=M(\ker V_F),
$$

where \(\operatorname{Stab}_{\mathrm{BP}}(F)\) denotes equality of the active
edge-length functions on the BP sign cover.

In this language, edge-sign separation says that the active edge-length
functions factor through the elementary-abelian character quotient recorded by
the labeled crossing characters \((e,C)\), and that this quotient is generic.
A generator that is mirror-compatible with \(e\) fixes \(\lambda_e\), while a
non-compatible crossing contributes a nonzero character coefficient generically.
Lemma 3 gives the local nonzero coefficient \(4h_ah_b\). The field-theoretic
task is to show that distinct labels give distinct quadratic characters in this
abelianized quotient of the BP cover.

This reformulation removes any dependence on an ordering of geometric
reflections. The labeled matrix records the abelian character quotient that is
visible to active edge lengths, even if the full iterated quadratic monodromy is
larger or nonabelian.

## Algebraic Theorem: Rank of the Projected Kernel

This part is unconditional.

Let

$$
\mathcal K_F=M(\ker V_F).
$$

Then

$$
\dim \mathcal K_F
=
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F).
$$

Proof.

Define

$$
T:A\to\mathbb F_2^{B_P}\oplus\mathbb F_2^{\mathcal L_F},
\qquad
T(\alpha)=(M\alpha,V_F\alpha).
$$

Let \(\pi\) be projection from \(\operatorname{im}T\) onto the violation
coordinates. Then

$$
\operatorname{im}\pi=\operatorname{im}V_F,
$$

so

$$
\operatorname{rank}\pi=\operatorname{rank}V_F.
$$

The kernel of \(\pi\) inside \(\operatorname{im}T\) is

$$
\ker\pi
=
\{(M\alpha,0)\mid V_F\alpha=0\}
=
\mathcal K_F\oplus\{0\}.
$$

By rank-nullity applied to
\(\pi:\operatorname{im}T\to\operatorname{im}V_F\),

$$
\dim\mathcal K_F
=
\dim\operatorname{im}T-\dim\operatorname{im}V_F
=
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F).
\quad\square
$$

## Theorem 2: Generic Labeled-Violation Rank Count

Assume \(K\ge2\), \(\Xi_F\neq\varnothing\), and the DDGP instance is generic in
the sense above. Then, outside a proper algebraic exceptional set,

$$
|\Xi_F|
=
2^{
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F)
}.
$$

Proof.

Let \(s_0\in\Xi_F\). By Theorem 1,

$$
\Xi_F=s_0\oplus\mathcal U_F(s_0).
$$

By Lemma 6,

$$
\mathcal U_F(s_0)=\mathcal K_F=M(\ker V_F).
$$

Therefore \(\Xi_F\) is an affine coset of the vector space
\(\mathcal K_F\), and

$$
|\Xi_F|
=
|\mathcal K_F|
=
2^{\dim\mathcal K_F}.
$$

Substituting the rank theorem gives the stated formula. \(\square\)

## Local Valid-Count Corollary

For a pruning-edge set \(P\), take

$$
F=E_0[L_P]\cup P.
$$

For \(K\ge2\), outside the same generic exceptional set, the number of local
valid DDGP solutions is

$$
|\Xi_P|
=
2^{
\operatorname{rank}
\begin{bmatrix}
M\\
V_F
\end{bmatrix}
-
\operatorname{rank}(V_F)
}.
$$

This is the form implemented by `scripts/ddgp_rank_count.py`.

## DMDGP Sanity Check

In the DMDGP, each predecessor clique is consecutive:

$$
U_i=\{v_{i-K},\dots,v_{i-1}\}.
$$

Dependency cones are suffixes of the order. The graph reflections become the
usual DMDGP partial reflections. A pruning edge kills exactly the suffix
reflections whose boundary separates the two endpoints. Therefore the labeled
matrix reduces to the SBBU symmetry test: a reflection survives exactly when no
active edge has an uncancelled crossing label.

Thus the DDGP rank formula contains the DMDGP symmetry count as the
consecutive-predecessor special case.

## Proof Status

Proved here:

- lateration choices are mirror reflections;
- mirror-compatible crossings have no single-edge mirror obstruction;
- a single non-compatible crossing changes edge length generically;
- cone masks span the local branch space;
- the labeled obstruction quotient \(Q_F\) makes the obstruction of a branch
  difference presentation-independent;
- outside a finite algebraic exceptional set, local solution codes are
  translates of the algebraic edge-preserving differences from a base code;
- primitive bracket obstructions are separated by a stretched-moment-curve
  specialization for \(K\ge2\);
- the generic labeled presentation lemma for \(K\ge2\);
- the rank identity for \(M(\ker V_F)\);
- the generic labeled-violation rank count theorem.

Remaining refinements:

- write the \(K=1\) degenerate case separately if it is ever needed;
- polish the block-reflection normal-form argument for publication, since the
  research note states it at the level needed for the current DDGP framework.

The active-edge formulation is stronger than the earlier base-filtered version:
it predicts both local base counts and local valid counts. In the current
generic stress tests it matches both counts in 2225/2225 stored cases for
DDGP/DMDGP instances in \(K=2\) and \(K=3\). An additional temporary \(K=2\),
\(n=13\), 80-instance multi-pruning stress run also matched both base and valid
counts in 80/80 cases.
