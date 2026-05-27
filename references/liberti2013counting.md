1

# Counting the number of solutions of $^{K}$ DMDGP instances

Leo Liberti$^{1,2}$, Carlile Lavor$^{3}$, Jorge Alencar$^{3}$, and Germano Abud$^{3}$

$^{1}$ IBM “T.J. Watson” Research Center, Yorktown Heights, 10598 NY, USA
leoliberti@us.ibm.com

$^{2}$ LIX, Ecole Polytechnique, 91128 Palaiseau, France
liberti@lix.polytechnique.fr

$^{3}$ Dept. of Applied Math. (IMECC-UNICAMP), State University of Campinas, 13081-970, Campinas – SP, Brazil
clavor@ime.unicamp.br, jorge.fa.lima@gmail.com, germano@famat.ufu.br

**Abstract.** We discuss a method for finding the number of realizations in $\mathbb{R}^K$ of certain simple undirected weighted graphs.

# 1 Introduction

In this paper we deal with Euclidean realizations of weighted graphs such that the Euclidean distance between pairs of realization points are the same as the weight on the corresponding edge. The vertex sets of our graphs are assumed to be ordered in certain ways formally described below. An order is given as a rank function $\rho$ mapping the vertex (having cardinality $n$) set into and onto the set $\bar{n} = \{1, \dots, n\}$. In general, for clarity of notation, we may identify a vertex with its rank, e.g. for any two vertices $u, v$ and an integer $K$, we write $u &lt; v$ or $v &gt; K$ to mean $\rho(u) &lt; \rho(v)$ and $\rho(v) &gt; K$.

The $K$-DISCRETIZABLE MOLECULAR DISTANCE GEOMETRY PROBLEM ($^{K}$ DMDGP) is as follows. Given a positive integer $K$, a simple undirected weighted graph $G = (V, E, d)$ where $d: E \to \mathbb{R}_+$, an order $&lt;$ on $V$ such that $\{u, v\} \in E$ for each $v &gt; K$ and $v - K \leq u \leq v - 1$, and a partial realization $\bar{x}: \{1, \dots, K\} \to \mathbb{R}^K$, does there exist a realization $x: V \to \mathbb{R}^K$ such that:

$$
\forall \{u, v\} \in E \quad \| x_u - x_v \| = d_{uv} \tag{1}
$$

and such that $x_v = \bar{x}_v$ for each $v \in \{1, \dots, K\}$? We remark that solving Eq. (1) for any given weighted graph and integer $K$ is known as the DISTANCE GEOMETRY PROBLEM (DGP). Both the DGP [18] and the $^{K}$ DMDGP [12] are NP-hard, even for fixed $K$.

The motivation for the name — molecular — stems from the natural application to finding molecular conformation (so $K = 3$). The vertices of the input graph $G$ are atoms, and the edges are pairs of atoms for which the distance is known. We focus on the important case of proteins: since all proteins consist of a backbone with some side chains, we consider the backbone as a natural vertex

---

order. Since covalent bond lengths are known, and the angles between covalent bonds is also known *[19]*, distances corresponding to pairs of atoms $\{v-1,v\}$ and $\{v-2,v\}$ are known. Nuclear Magnetic Resonance (NMR) experiments provide an estimation of distances shorter than around 6Å, which covers the case of pairs $\{v-3,v\}$ (as well as other pairs — the backbone folds in space, and it often happens that two atoms that are far apart in the order are actually close in Euclidean space) *[19]*. Moreover, for elementary geometrical reasons, it is always possible to fix positions of the first, second and third atom in the protein backbone so that the inter-atomic distances over $\{1,2,3\}$ are satisfied. Thus, protein backbones provide natural examples of ^{∗}DMDGP instances *[9]*. In the following, we shall partition the edge set $E$ into the discretization edges $E_{D}=\{\{v,v-j\}\mid j\in\{1,\ldots,K\}\}$ and the pruning edges $E_{P}=E\smallsetminus E_{D}$. We let $m=|E|$.

Finding realizations for general graphs usually involves a continuous search *[13]*, but if the graph is rigid *[6]* then a discrete search type is possible *[10]*. It was observed in *[8]* that ^{∗}DMDGP graphs are Henneberg graphs, which are known to be rigid *[20]*. In *[11]* we proposed a discrete search algorithm called Branch-and-Prune (BP), where the discretization edges are used to make sure that only a discrete set of points needs to be checked for feasibility w.r.t. Eq. (1), and the pruning edges are used to reduce the search space.

Since every vertex $v>K$ is adjacent to (at least) its $K$ immediate predecessors, if we know the position $x_{u}$ of each of these predecessors $u$ of $v$, then $x_{v}$ is at the intersection of $K$ spheres in $\mathbb{R}^{K}$ *[2]*. Provided the strict simplex inequalities (a generalization of the triangle inequalities to $\mathbb{R}^{K}$ *[8]*) hold, this intersection is either empty or consists of exactly two points. This provides an inductive step to find the next vertex in the order, having placed all its predecessors. The base case is dealt with since we are given the partial realization $\bar{x}$. Since at each step there may be two feasible positions $x_{v}$ for the next vertex $v$, in the worst case the BP yields an exponentially large search tree, where each node $x_{v}$ at level $v$ is a possible position for vertex $v$. Since the first branch occurs at level $K+1$, this worst-case tree has $2^{n-K}$ leaf nodes. Each leaf node $x_{n}$ corresponds to a unique path from the root $x_{1}$ to $x_{n}$, which therefore encodes a valid realization $x=(x_{1},\ldots,x_{n})$. We let $X$ be the set of all these realizations.

In this paper, we propose an efficient method for computing the cardinality of $X$.

## 2 Motivation

Knowing $|X|$ is important for at least two practical reasons. First, in the application of DGP to proteomics, the set $X$ is of interest to biochemists, who will evaluate each potential backbone according to chemical criteria. If $|X|$ is too large, this evaluation might be too costly; on the other hand, if $|X|$ is too small, $X$ might not contain the “correct” backbone. This observation might sound strange to mathematicians, but one must not forget that the ^{∗}DMDGP provides a model of reality, rather than being reality itself: none of the realizations in $X$

---

might be really correct from the point of view of the biochemical practitioner, but some may be close enough for him or her to recognize them. From a different point of view, the experimental data set usually contains errors *[1]* which might influence the number of realizations in $X$: a small $|X|$ might be evidence of wrong data.

Secondly, the class of globally rigid (also known as “uniquely realizable” *[7]*) graphs, i.e. those for which it can be shown that $|X|=1$, is interesting because in several DGP applications, such as e.g. to wireless sensor localization, ensuring that sufficient distance data are known for the the graph having a unique realization is of paramount importance: recovering a large set of different possible networks obviously prevents practitioners from understanding the actual network geometry. Necessary (combinatorial) conditions for a graph to be globally rigid are given in *[7]*: informally speaking, if removing a certain edge from a rigid graph still yields a rigid graph, the edge is redundant; in a redundantly rigid graph, all edges are redundant. Redundant rigidity turns out to be a necessary condition for unique realizability. Although there are exact methods for verifying whether a graph is redundantly rigid for $K\in\{1,2\}$, no such method is known for higher dimensions. A randomized $O(n^{2}m)$ method is given in *[7]*.

Although no necessary and sufficient conditions for unique realizability is known so far, several different sufficient conditions are known. Cliques are obviously globally rigid, and the realization can be found in polynomial time *[3]*. Trilateration graphs are those for which there exists a vertex order where each $v>K$ has at least $K+1$ adjacent predecessors: these can be shown to have a unique realization, which can be found in polynomial time *[4]*. The graphs occurring in ^{∗}DMDGP are a natural generalization of trilateration graphs, insofar as they require at least $K$ adjacent predecessors. As shown in *[14]*, in general such graphs are not globally rigid, but the number of realizations can be counted in time $O(n+m)$, as shown in Sect. 4 below; so those ^{∗}DMDGP graphs that are globally rigid can be recognized in polynomial time (under some genericity assumptions, see Sect. 3.2). Uniquely localizable graphs possess a unique realization in a given $K$, and no other realization for any higher value of $K$. It is shown in *[15]* that these graphs can be realized in polynomial time (up to some approximation constant) by solving a semidefinite programming problem.

## 3 Background material

Although our method for computing $|X|$ is straightforward, is rests on many known but nontrivial results, which we summarize here.

### 3.1 Incongruence

Two sets of points in $\mathbb{R}^{K}$ are congruent if there is a sequence of translations, rotations and reflections that turns one into the other. Since any realizable graph has uncountably many congruent realizations, we are only interested in the number of incongruent ones. Unfortunately, the way we defined $X$ above (i.e. $X$ is

---

the set of solutions found by the BP algorithm on ^{e}DMDGP instances) is only partially correct in this respect. Because the realizations of Henneberg graphs are rigid frameworks, each realization in $X$ is rigid; so the fact that the first $K$ vertices are fixed in given positions $\bar{x}_{1},\ldots,\bar{x}_{K}$ eliminates rotations and translations. By Thm. *[9, Thm. 2]*, with $K=3$ there is a “fourth-level symmetry” in $X$: half of the realizations in $X$ are reflections of the other half along the plane through $\bar{x}_{1},\bar{x}_{2},\bar{x}_{3}$. This was generalized in *[14]* for any $K$.

So that the definition of $X$ is consistent with $X$ being a set of incongruent realizations, we simply modify the BP algorithm to choose any of the two possible positions for $x_{K+1}$ (without exploring the other), and start branching from level $K+2$.

### 3.2 Probability 1

The theory supporting the BP algorithm is always based on the edge weight function $d$ satisfying the strict simplex inequalities (i.e. the Cayley-Menger determinant of each $K$-subsequence of vertices in the given order multiplied by $(-1)^{K+1}$ is strictly positive). Otherwise, the intersection of $K$ spheres in $\mathbb{R}^{K}$ might have uncountable cardinality, or be a singleton set. These occurrences only happen when the ^{e}DMDGP instance is YES, and the values assigned to $d$ yield zero Cayley-Menger determinants *[8]*, i.e. they satisfy a certain given system of polynomial equations. Such systems define manifolds of Lebesgue measure zero in $\mathbb{R}^{K}$. Moreover, it is easy to prove that all points in all realizations in $X$ are in a ball centered at $\bar{x}_{1}$ with radius bounded by the sum of all edge distances. So, the probability of uniformly sampling $d$ satisfying these equations is zero. This in turn means that the probability of uniformly sampling $d$ such that it yields a YES ^{e}DMDGP instance satisfying the strict simplex inequalities is 1. Accordingly, we state most of our results “with probability 1”.

There are at least three related concepts in the literature. The first, genericity (in the standard sense), requires that there should be no rational polynomial satisfying the instance data $d$. This condition is “too strong”, in the sense that it would require at least one value of $d$ to be transcendental, which makes little sense for computers. The second concept requires that all minors of the complete rigidity matrix are nontrivial *[5]*. The third requires that $d$ is a rational function contained in the (open) complement of the set of those rational functions $d^{\prime}$ yielding zero Cayley-Menger determinants *[16, 17]*. The notion we employ is very similar to both the second and the third concept.

### 3.3 Partial reflections

For any realization $x\in X$ and $v\in V$ with $v>K$, we let $R^{v}_{x}$ be the reflection along the hyperplane through $x_{v-K},\ldots,x_{v-1}$, as shown in Fig. 1. Now, for any $v>K$, we define a partial reflection operator with respect to $x$ as:

$g_{v}(x)=(x_{1},\ldots,x_{v-1},R^{v}_{x}(x_{v}),R^{v}_{x}(x_{v+1}),\ldots,R^{v}_{x}(x_{n})).$ (2)

---

![img-0.jpeg](img-0.jpeg)
Fig. 1. The action of the reflection  $R_{s}^{v}$  in  $\mathbb{R}^K$ .

The partial reflection  $g_v$  acts on a realization  $x$  by reflecting all vectors from rank  $v$  onwards. We define a product between partial reflections by setting  $g_u g_v = g_u \circ g_v$  for all  $u, v &gt; K$ , i.e.  $g_u g_v$  is the operation consisting in applying  $g_v$  first, and then  $g_u$  later to a realization  $x \in X$ . More precisely, for  $v &gt; u &gt; K$  and  $x \in X$ ,

$$
\begin{array}{l} g _ {u} g _ {v} (x) = g _ {u} \left(g _ {v} (x)\right) = g _ {u} \left(x _ {1}, \dots , x _ {v - 1}, R _ {s} ^ {v} \left(x _ {v}\right), \dots , R _ {s} ^ {v} \left(x _ {n}\right)\right) = \\ = \left(x _ {1}, \dots , x _ {u - 1}, R _ {s} ^ {u} \left(x _ {u}\right), \dots , R _ {s} ^ {u} \left(x _ {v - 1}\right), R _ {g _ {v} (x)} ^ {u} \left(x _ {v}\right), \dots , R _ {g _ {v} (x)} ^ {u} \left(x _ {n}\right)\right) \\ \end{array}
$$

(the case for  $u &lt; v$  is similar). Notice that the action of the left operand  $g_{u}$  after rank  $v$  does not apply  $R_{s}^{u}$  to the components of the argument, but  $R_{g_{v}(x)}^{u}$ . By [12, Lemma 2], this product is commutative.

Now let  $\Gamma_{D} = \{g_{v} \mid v &gt; K\}$ , and consider the set  $\mathcal{G}_D = \langle \Gamma_D \rangle$  generated by all possible products of elements in  $\Gamma_{D}$ . By [12],  $\mathcal{G}_D$  turns out to be the invariant group of the set of realizations  $X_{D}$  consisting of all the possible realizations found by the BP algorithm on the graph  $G_{D} = (V, E_{D})$  induced by the discretization edges (see Fig. 2). Our purpose is to find the invariant group  $\mathcal{G}_P$  of the set of realizations of the given graph  $G$ , which we assume to have a nontrivial set of  $E_{P}$  of pruning edges. Let the span of a pruning edge  $\{u, w\} \in E_P$  be the set  $S^{uw} = \{u + K + 1, \dots, w\}$  (assuming  $u &lt; w$ ; if  $w &gt; u$  we let  $S^{uw} = S^{wu}$ ). By [12],  $\mathcal{G}_P$  is the subgroup of  $\mathcal{G}_D$  generated by

$$
\Gamma_ {P} = \left\{g _ {v} \mid v &gt; K \wedge \forall \{u, w \} \in E _ {P} (v \notin S ^ {u w}) \right\}. \tag {3}
$$

In other words, only those vertices that are not in the span of any pruning edge give rise to partial reflection operators that generate the discretization group  $\mathcal{G}_P$ .

---

![img-1.jpeg](img-1.jpeg)
Fig. 2. On the left: the set  $X_{D}$  of realizations of the graph induced by the discretization edges. On the right: the effect of the pruning edge  $\{1,4\}$  on  $X_{D}$ .

# 4 Counting incongruent realizations

By [12, Thm. 4], there is an integer  $\ell$  such that  $|X| = 2^{\ell}$  with probability 1. We can easily refine the proof of this result so that it says something more precise on  $\ell$ .

Proposition 1. With probability 1,  $|X| = 2^{|\Gamma_P|}$ .

Proof. The following statements hold with probability 1. By [12],  $\mathcal{G}_D \cong C_2^{n - K}$  (where  $C_2$  is the cyclic group of order 2), so that  $|\mathcal{G}_D| = 2^{n - K}$ . Since  $\mathcal{G}_P \leq \mathcal{G}_D$ ,  $|\mathcal{G}_P|$  divides the order of  $|\mathcal{G}_D|$ . By elementary group theory,  $|\mathcal{G}_P| = 2^{|\Gamma_P|}$ . By Thm. [14, Thm. 6.4], the action of  $\mathcal{G}_P$  on  $X$  only has one orbit, i.e.  $\mathcal{G}_Px = X$  for any  $x \in X$ . We remark that every partial reflection operator is idempotent, i.e.  $g_v^2 = 1$ , and hence  $g_v^{-1} = g_v$  for all  $v &gt; K$ . Thus, if  $gx = g'x$  for two  $g, g' \in \mathcal{G}_P$  and  $x \in X$ , then  $(g')^{-1}gx = x$ , which implies  $g'gx = x$ , which implies  $g'g = 1$  whence  $g' = g$ . This means that  $|\mathcal{G}_Px| = |\mathcal{G}_P|$ . Thus, for any  $x \in X$ ,  $|X| = |\mathcal{G}_Px| = |\mathcal{G}_P| = 2^{|\Gamma_P|}$ .

Now, all that remains to do is to present an algorithm to compute  $|\Gamma_P|$ . This follows directly from the definition in Eq. (3). We let  $b = (b_{K + 1},\ldots ,b_n)$  be an array initialized so that  $b_{i} = 1$  for all  $i$  in  $\{K + 1,\dots ,n\}$ . Then we scan every edge  $\{u,v\}$  in  $E_P$ , and for each  $i$  in  $S^{uv}$  we set  $b_{i} = 0$ . Finally,  $|\Gamma_P| = \sum_{i = K + 1}^n b_i$ . This algorithm runs in  $O(n + m)$ . We remark that, by Sect. 3.1, if  $X$  is required to only contain incongruent realizations, then the first component of  $b$  should be  $b_{K + 2}$  rather than  $b_{K + 1}$ .

Obviously, if  $|\Gamma_P| = 1$ , then the *DMDGP graph is globally rigid (with probability 1).

# Acknowledgments

Financial support is gratefully acknowledged from French National Research Agency project ANR-10-BINF-03-08 "Bip:Bip", and the Brazilian research agencies FAPESP, CNPq and CAPES.

---

References

- [1] Berger, B., Kleinberg, J., Leighton, T.: Reconstructing a three-dimensional model with arbitrary errors. Journal of the ACM 46(2) (1999) 212–235
- [2] Coope, I.: Reliable computation of the points of intersection of $n$ spheres in $\mathbb{R}^{n}$. Australian and New Zealand Industrial and Applied Mathematics Journal 42 (2000) C461–C477
- [3] Dong, Q., Wu, Z.: A linear-time algorithm for solving the molecular distance geometry problem with exact inter-atomic distances. Journal of Global Optimization 22 (2002) 365–375
- [4] Eren, T., Goldenberg, D., Whiteley, W., Yang, Y., Morse, A., Anderson, B., Belhumeur, P.: Rigidity, computation, and randomization in network localization. IEEE Infocom Proceedings (2004) 2673–2684
- [5] Graver, J.: Rigidity matroids. SIAM Journal on Discrete Mathematics 4 (1991) 355–368
- [6] Graver, J., Servatius, B., Servatius, H.: Combinatorial Rigidity. American Mathematical Society (1993)
- [7] Hendrickson, B.: Conditions for unique graph realizations. SIAM Journal on Computing 21(1) (1992) 65–84
- [8] Lavor, C., Lee, J., John, A.L.S., Liberti, L., Mucherino, A., Sviridenko, M.: Discretization orders for distance geometry problems. Optimization Letters 6 (2012) 783–796
- [9] Lavor, C., Liberti, L., Maculan, N., Mucherino, A.: The discretizable molecular distance geometry problem. Computational Optimization and Applications 52 (2012) 115–146
- [10] Lavor, C., Liberti, L., Maculan, N., Mucherino, A.: Recent advances on the discretizable molecular distance geometry problem. European Journal of Operational Research 219 (2012) 698–706
- [11] Liberti, L., Lavor, C., Maculan, N.: A branch-and-prune algorithm for the molecular distance geometry problem. International Transactions in Operational Research 15 (2008) 1–17
- [12] Liberti, L., Lavor, C., Mucherino, A.: The discretizable molecular distance geometry problem is easier on proteins. In Mucherino, A., Lavor, C., Liberti, L., Maculan, N., eds.: Distance Geometry: Theory, Methods, and Applications. Springer, New York (2013)
- [13] Liberti, L., Lavor, C., Mucherino, A., Maculan, N.: Molecular distance geometry methods: from continuous to discrete. International Transactions in Operational Research 18 (2010) 33–51
- [14] Liberti, L., Masson, B., Lavor, C., Lee, J., Mucherino, A.: On the number of realizations of certain Henneberg graphs arising in protein conformation. Discrete Applied Mathematics (accepted)
- [15] Man-Cho So, A., Ye, Y.: Theory of semidefinite programming for sensor network localization. Mathematical Programming B 109 (2007) 367–384
- [16] Nie, J., Ranestad, K., Sturmfels, B.: The algebraic degree of semidefinite programming. Mathematical Programming A 122(2) (2010) 379–405
- [17] Ranestad, K., Sturmfels, B.: Personal communication (2013)
- [18] Saxe, J.: Embeddability of weighted graphs in $k$-space is strongly NP-hard. Proceedings of 17th Allerton Conference in Communications, Control and Computing (1979) 480–489

---

19. Schlick, T.: Molecular modelling and simulation: an interdisciplinary guide. Springer, New York (2002)
20. Tay, T.S., Whiteley, W.: Generating isostatic frameworks. Structural Topology 11 (1985) 21–69

View publication stats