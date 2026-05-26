Discrete Applied Mathematics 165 (2014) 213-232

ELSEVIER

Contents lists available at ScienceDirect

Discrete Applied Mathematics

journal homepage: www.elsevier.com/locate/dam

SCIENCE DIRECT MATHEMATICS

# On the number of realizations of certain Henneberg graphs arising in protein conformation*

Leo Liberti $^{a,b,*}$, Benoit Masson $^{c}$, Jon Lee $^{d}$, Carlile Lavor $^{e}$, Antonio Mucherino $^{c}$

$^{a}$ LIX, École Polytechnique, 91128 Palaiseau, France
$^{b}$ IBM "T.J. Watson" Research Center, New York, USA
$^{c}$ IRISA, University of Rennes 1, Rennes, France
$^{d}$ IOE, University of Michigan, Ann Arbor (MI), USA
$^{e}$ Department of Applied Mathematics (IMECC-UNICAMP), University of Campinas, 13081-970, Campinas - SP, Brazil

# ARTICLE INFO

Article history:

Received 16 September 2011

Received in revised form 8 November 2012

Accepted 24 January 2013

Available online 6 March 2013

Keywords:

Distance geometry

Graph rigidity

Branch-and-Prune

Partial reflection

Protein conformation

# ABSTRACT

Several application fields require finding Euclidean coordinates consistent with a set of distances. More precisely, given a simple undirected edge-weighted graph, we wish to find a realization in a Euclidean space so that adjacent vertices are placed at a distance which is equal to the corresponding edge weight. Realizations of a graph can be either flexible or rigid. In certain cases, rigidity can be seen as a property of the graph rather than the realization. In the last decade, several advances have been made in graph rigidity, but most of these, for applicative reasons, focus on graphs having a unique realization. In this paper we consider a particular type of weighted Henneberg graphs that model protein backbones and show that almost all of them give rise to sets of incongruent realizations whose cardinality is a power of two.

© 2013 Elsevier B.V. All rights reserved.

# 1. Introduction

The fundamental problem of Distance Geometry (DG) is that of determining the Euclidean coordinates corresponding to a given set of distances [9]. By "distances" we mean here a set $E$ of unordered pairs $\{u, v\}$ of vertices of a set $V$ together with a function $d: E \to \mathbb{R}_+$ mapping each pair $\{u, v\}$ to a distance value $d_{uv}$ between $u$ and $v$. We thus formalize this problem as the

DISTANCE GEOMETRY PROBLEM (DGP). Given a positive integer $K$ and a weighted simple undirected graph $G = (V, E, d)$, where $d: E \to \mathbb{R}_+$, find a function $x: V \to \mathbb{R}^K$ such that

$$
\forall \{u, v \} \in E \quad \| x _ {u} - x _ {v} \| _ {2} = d _ {u v}. \tag {1}
$$

A vertex function $x: V \to \mathbb{R}^K$ is a realization. A realization is valid if it satisfies Eq. (1). The DGP, also known as the EMBEDDABILITY problem, is strongly NP-hard even when $K = 1$ and the edge weights are restricted to take values in

0166-218X/$ - see front matter © 2013 Elsevier B.V. All rights reserved.

doi:10.1016/j.dam.2013.01.020

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

![img-0.jpeg](img-0.jpeg)
Fig. 1. A $^{K}$ DMDGP graph with $K = 3$.

$\{1,2\}$ [42]. It was shown in [1] that all graphs having a realization in some Euclidean space are realizable when $K \geq \frac{1}{2} (\sqrt{8|E| + 1} - 1)$, so, in a sense, the smaller the dimension, the harder the problem.

Historically, DG was first studied by Menger in the 1930s [31] and then fully discussed and extended by Blumenthal [4]. At that time, the main notion was that of a distance space $(S,d)$, where $S$ is a set and $d:S\times S\to \mathbb{R}$ is a distance function. The "fundamental problem" (also called subset problem by Blumenthal) was that of determining necessary and sufficient conditions for given weighted sets $(S,d)$ to be distance spaces. In other words, given a square $n\times n$ matrix $D = (d_{ij})$, determine whether it is a distance matrix, i.e. a matrix for which there exists a realization $x$ in $K$ dimensions such that $d_{ij} = \| x_i - x_j\|$ for all $i,j\in S$. Schoenberg [44] proved in 1935 that this is the case if and only if the matrix $-V^{\top}DV$ is Positive Semi-Definite (PSD), where $V$ is the $n\times (n - 1)$ matrix:

$$
\frac {1}{\sqrt {2}} \left( \begin{array}{c c c c} - 1 &amp; - 1 &amp; \dots &amp; - 1 \\ 1 &amp; 0 &amp; \dots &amp; 0 \\ 0 &amp; 1 &amp; \dots &amp; 0 \\ \vdots &amp; \vdots &amp; \ddots &amp; \vdots \\ 0 &amp; 0 &amp; \dots &amp; 1 \end{array} \right).
$$

We remark that Schoenberg used another notation; our notation is taken from [10, Section 5.4.1]. Since a matrix can be ascertained to be PSD in polynomial time, Blumenthal's subset problem is in $\mathbf{P}$.

A finite simple weighted undirected graph may have uncountably many or finitely many realizations. In this paper we are interested in rigid graphs, which have finitely many realizations. Graph rigidity was discussed time and again in connection with several different application fields, the main ones being statics, molecular biology, robotics and localization of wireless sensor networks. Historically, this generated several slightly different definitions of the concept of graph rigidity. An advanced but didactical account can be found in [14].

We consider the case where the graph $G$ has a vertex order such that each $(K + 1)$-tuplet of consecutive vertices of $V$ induces a $(K + 1)$-clique as a subgraph of $G$ (see Fig. 1). The class of DGP instances with such an order, collectively known as the generalized DISCRETIZABLE MOLECULAR DISTANCE GEOMETRY PROBLEM (K DMDGP), is also NP-hard when $K = 3$ [22]. The K DMDGP order is a particular type of Henneberg order [46], and therefore characterizes a rigid graph. However, in most cases the graph is not uniquely rigid. In this paper we prove that, for almost all edge weight functions, K DMDGP graphs have a number of incongruent realizations which is a power of two. As Hendrickson pointed out in [16, Section 3], non-uniquely rigid graphs are not well studied. This might be because one of the most important applications of DG (the localization of sensor networks [12]) requires the input data to be dense enough so that the graph is uniquely rigid [3]. More importantly, Hendrickson also pointed out that the key to non-unique rigidity is partial reflection of the realization: we shall make use of this concept to derive our result.

Our motivation for studying rigid graphs which are not uniquely rigid arises from the conformation of proteins. In general, a molecule can be seen as a unit sphere graph [19] defined by a distance threshold given by the resolution scope of Nuclear Magnetic Resonance (NMR) machinery (this is between $5\AA$ and $6\AA$ [43]). In other words, all atom pairs closer than this threshold are adjacent. The class of DGP instances defined by these unit sphere graphs, with $K = 3$, is known as the MOLECULAR DISTANCE GEOMETRY PROBLEM (MDGP) (see [29] and citations therein). We argued in [28,22] that whenever the instance represents a protein backbone, the natural backbone order is a $^{K}$ DMDGP order. Since, however, the DGP and its variants require exact distances, whereas NMR data are best represented by intervals due to their associated measurement error [2], a remark is in order. Certain interatomic quantities are known reasonably precisely: notably, covalent bonds and covalent angles [43]. Moreover, because of the scaling between interval width and NMR machinery resolution scope, certain NMR intervals can be considered as finite sets of values [39]. This allows one to define a virtual order on protein backbones such that each consecutive quadruplet of atoms induces a 4-clique with the property that at most one of its edges is weighted by a finite set of values (instead of a single value) [23]. Under these conditions, the resulting graph can be seen as a collection of finitely many $^{K}$ DMDGP graphs, each of which is rigid.

The present work falls into a sequence of works about solving realizability problems using the BP approach. In this paper we make an important theoretical contribution about the structure of the solution set of $^{K}$ DMDGP instances: given any

---

solution x, the others can all be obtained as “partial reflection” operators applied to x. This can also be used to speed up the time to find X, as shown in [33], where we apply some of the ideas of this paper to the problem of finding 3D realizations of protein backbones.

The rest of this paper is organized as follows. We recap some rigidity definitions and results (Section 2) and sketch the Branch-and-Prune (BP) algorithm for finding^{K}DMDGP graph realizations (Section 3). We then give formal descriptions of the^{K}DMDGP and of the BP algorithm (Section 4). Next, we study the geometrical properties of the BP search tree (Section 5), and prove that the number of solutions of YES instances of the^{K}DMDGP is a power of two with probability one (Section 6). We exhibit a (zero measure) family of counterexamples to the “power of two” conjecture in Section 7, and finally we extend our results to the application-specific setting of protein conformation (Section 8).

## Rigidity

Two distance spaces (S, d) and (S′, d′) in ℝ^{K} are congruent when there exists an affine operator *T*:ℝ^{*K*}→ℝ^{*K*} such that the restriction of T to S is a bijection *S*→*S*^{′} and *d*(*p*,*q*) = *d*^{′}(*T**p*,*T**q*) for all *p*, *q* ∈ *S*. Such operators are also called isometries. A framework is a pair (G, x) where *G* = (*V*, *E*, *d*) is a simple undirected graph weighted by *d*:*E*→ℝ_{+} and *x*:*V*→ℝ^{*K*} is a valid realization of G. A displacement of a framework (G, x) is a continuous function *y*:[0, 1]→ℝ^{*K*|*V*^{′}} such that *y*(0) = *x* and *y*(*t*) is a valid realization of G for all *t* ∈ [0, 1]. A flexing of a framework (G, x) is a displacement y of (G, x) such that *y*(*t*) is incongruent to x for any *t* ∈ (0, 1]. A framework is flexible if it has a flexing, otherwise it is rigid. The rigidity matrix R of a framework (G, x) has |E| rows and K|V| columns; the row indexed by {*u*,*v*} ∈ *E* has exactly 2K nonzero components, namely *x*_{*u**k*} - *x*_{*v**k*} in the columns indexed by (*u*,*k*) and *x*_{*v**k*} - *x*_{*u**k*} in the columns indexed by (*v*,*k*) (for *k* ≤ *K*). The complete rigidity matrix R refers to the case where G is the full clique on V. The rigidity matrix is used to give a stronger definition of rigidity (namely that of infinitesimal rigidity) which is used in statics [47,48]. Intuitively, infinitesimally rigid structures are resilient to collapse when certain forces are applied to them. Since for the kinds of frameworks considered in this paper rigidity and infinitesimal rigidity are equivalent [46, p. 23], we only limit the discussion to rigidity.

As defined above, rigidity is a property of frameworks rather than of graphs. It turns out, however, that if a graph has a certain type of rigid framework, then all its frameworks are rigid. More precisely, if x is a valid realization whose components are all algebraically independent over ℚ (i.e. there is no polynomial over ℚ having all the components as roots), then x is called a generic realization: if (G, x) is rigid, then all generic realizations of G are rigid, and we can then say that G is a rigid graph. This statement implicitly assumes that G is not weighted: in other words, a graph is rigid if its rigidity properties depend on the graph topology itself rather than the edge weight function, apart perhaps from a set of edge weight functions which have at most the cardinality of algebraic relations over ℚ. Since the latter is a countable cardinality, whereas the set of all edge weight functions is uncountable, that is another way of saying that almost all edge weight functions yield a rigid framework if the graph is rigid. Having said that, the definition of generic realization given above is too stringent. Typically, rigid graphs might fail to yield rigid frameworks for edge weight functions which cause (*K* + 1)-cliques to be embedded in a space of dimension < *K* (e.g. a 3-cycle embedded in ℝ^{2} as three collinear points). Graver's definition of a generic realization [13] makes this concept precise: each nontrivial minor of the complete rigidity matrix must be nonzero. This implies that there is no need for calling in the very strong requirement of full algebraic independence over ℚ: it suffices to make sure that the components of the realization giving rise to a rigid framework do not satisfy the polynomial equations engendered by the nontrivial minors of the complete rigidity matrix. If G has such a realization, then the graph is rigid. Also see Section 3.1.

Graphs with a unique realization are known as uniquely or globally rigid. Global rigidity recently generated several theoretical advances, such as a polynomial algorithm for realizing such graphs using the duality theory of Semi-Definite Programming (SDP) [45] and the equivalence between unique K-localizability and generic universal rigidity [49] (among others). Unique K-localizability and universal rigidity are strong forms of unique realizability, which also involve uniqueness in higher dimensions. All these interesting properties arise in the study of localization of sensor networks, and are based on two application-driven requirements: unique realizability and the presence of a set of vertices (called beacons or anchors) whose realization is known a priori. Since at any time the network has exactly one localization to be determined, graphs should be dense enough so that they guarantee some unique realizability property. Moreover, most sensor networks are linked to other communication networks by means of routers (representing the anchors) whose position in space is fixed and known. When there are sufficiently many anchors in general positions, one can guarantee that any valid realization of the network can only take place in the full K space, and not in any smaller dimensional one. The application of rigidity to protein fails to provide either anchors or unique realizability guarantees. Quite on the contrary, it is interesting, from a biological point of view, to have a list of possible bio-polymers that a given set of distances can realize.

### Henneberg type I graphs

Several important results on the rigidity of frameworks date from the end of the XIX century [8,41], motivated by the construction of buildings whose supporting structures consisted of bars and joints. It appears from the literature of the period that verifying rigidity was an inductive process: one would start with a rigid structure and then add rigid components to it, so that the resulting structure would also be rigid. Henneberg [17] was the first to formalize two inductive steps for verifying

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

![img-1.jpeg](img-1.jpeg)
Fig. 2. General case for the intersection $P$ of three spheres in $\mathbb{R}^3$.

rigidity. The first of these (known as Henneberg type I step [46,18]) can be paraphrased (and generalized) as follows: if there is an order on $V$ such that the first $K$ vertices have a known realization, and such that every subsequent vertex is adjacent to at least $K$ predecessors, then the graph almost certainly has a rigid realization in $\mathbb{R}^K$. This idea was already present in the works of Saviotti [40, Section XI, p. 57 and Fig. 30, pl.XV]. A similar order, which asks for the first $K$ vertices to be a clique and each subsequent vertex to be adjacent to at least $K$ predecessors, is called discretization vertex order, and the problem of finding such an order given a graph is known as the DISCRETIZATION VERTEX ORDER PROBLEM (DVOP) [20]. The DVOP is NP-complete, but only because the $K$-clique problem trivially reduces to it. Once the initial $K$-clique is known, a greedy procedure can find the order or decide it does not exist. The DVOP is therefore in $\mathbf{P}$ for fixed $K$. Because in Henneberg type I orders a partial realization is given for the first $K$ vertices, it follows that the initial $K$-clique is known in advance; finding such an order in a graph is therefore also in $\mathbf{P}$. Graphs with a Henneberg type I vertex order are also called Henneberg type I graphs. It is shown in [34] that graphs with a discretization vertex order are rigid; the realization problem for such graphs is known as the DISCRETIZABLE DISTANCE GEOMETRY PROBLEM (DDGP) and is NP-hard [34]. We remark that requiring that each vertex should have at least $K + 1$ adjacent predecessors (instead of only $K$) yields a $K$-trilateration order: graphs with such orders are called $K$-trilateration graphs and can be realized in polynomial time [12]. It is interesting that the difference between the definitions of Henneberg type I orders and $K$-trilateration orders is as small as possible, and yet marks the distinction between an easy and a hard corresponding realization problem.

In the rest of this paper we are concerned with a particular type of Henneberg type I orders, namely those that ensure that the next vertex is adjacent to exactly $K$ immediate predecessors, and perhaps also to other (not necessarily immediate) ones. We call these $^K$ DMDGP orders: they are important in using distances to find the conformation of proteins in space. As mentioned above, this restriction on the order does not make the corresponding realization problem any easier from the worst-case computational complexity standpoint [22]. It does guarantee rigidity, however, and it allows us to devise a recursive procedure for finding all incongruent realizations. Moreover, as we shall see, it also allows us to determine that the number of these realizations is a power of two for almost all edge weight functions.

# 3. Sphere intersections and Branch-and-Prune

Henneberg type I graphs are rigid because they enjoy the Sphere Intersection Property (SIP). Let $v \in V$ and suppose the first $v - 1$ vertices have already been realized in $\mathbb{R}^K$: we know then that the $v$-th vertex is adjacent to at least $K$ predecessors. Consider $K$ of these. Since they precede $v$ in the order, their position in $\mathbb{R}^K$ is already known by induction. Their distance to $v$ is known because they are adjacent to $v$. Therefore, the position for $v$ can only belong to the intersection of $K$ spheres centered at these $K$ adjacent predecessors. Since the intersection of $K$ spheres in $\mathbb{R}^K$ generally consists of either 0 or 2 points (see Fig. 2), no vertex can ever be placed in an uncountable set of positions in $\mathbb{R}^K$, as the definition of flexing would imply. Here, the term generally has a meaning similar to the one given by Graver's definition of generic rigidity [13]: the set of edge weight functions for which the statement fails to hold has countable cardinality. Since such sets have Lebesgue measure 0 in $\mathbb{R}^K$, we also speak of this statement holding with probability 1. In this particular case, the set of $K$-tuples of spheres for which the SIP does not hold has Lebesgue measure 0 in the set of all possible $K$-tuples of spheres (see Fig. 3). See Section 3.1 for a more complete discussion about probability 1 statements.

The SIP can be exploited algorithmically by a recursive algorithm which will test in turn each of the two possible positions for the current vertex. If they yield valid partial realizations, then the algorithm calls itself recursively and tries to place the next vertex. This algorithm is called Branch-and-Prune (BP) [28]. The branching occurs because at each rank in the vertex order, two positions are possible. Pruning occurs when a vertex is adjacent to more than $K$ adjacent predecessors: when a proposed position for the $v$-th vertex is inconsistent with some adjacent predecessor, recursion on the next vertex does not occur. The BP algorithm was originally only defined for $^K$ DMDGP graphs [21,22], but was then extended to several other settings: for Henneberg type I graphs [34], for certain types of interval-weighted graphs related to proteins [32,26,23], and for the purpose of overcoming a technical limitation of NMR machinery, which generally only provides reliable distance measures for pairs of hydrogen atoms [35,24,25,27]. A publicly available BP implementation is described in [38]. The

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

![img-2.jpeg](img-2.jpeg)
Fig. 3. When three adjacent predecessors are collinear in  $\mathbb{R}^3$ , the sphere intersection (the thick black circle) might have uncountable cardinality.

current computational state-of-the-art for the BP algorithm is attained with a parallel BP implementation [36,37], which can realize a protein backbone of  $10^{4}$  atoms (183 444 distances) in  $\mathbb{R}^3$  in 1.57 s of CPU time on a cluster of 64 nodes.

It is easy to see that the BP yields a worst-case exponential behavior, occurring when each vertex has exactly  $K$  adjacent predecessors. In such a case, the BP search tree is a full binary tree of height at most  $|V|$  and width  $2^{|V| - K}$  attained at the last level. Paths of length  $|V|$  from the root to a leaf node encode realizations of the input graph and hence denote YES instances of the  $^{\mathrm{K}}$ DMDGP, whereas a tree with height strictly less than  $|V|$  certifies a NO instance. In practice, however, the BP outperforms its continuous search competitors in both efficiency and reliability [22]. One particularly useful feature of BP is that, because the search is complete, it finds the set  $X$  of all incongruent realizations for a given graph, whereas most other DGP algorithms only find one realization. As already remarked, this is useful in biology because it allows one to list all the bio-polymers that are consistent with a particular set of distances.

In all our computational tests on  $^{\mathrm{K}}$ DMDGP instances, we observed that the number of incongruent realizations is a power of two: this comes to no surprise in the exponential worst case mentioned above, but there is no apparent reason why this should be the case when adjacent predecessors also include other vertices than the  $K$  immediate predecessors; and, indeed, in Section 7 we exhibit a set of counterexamples to the conjecture that all YES instances of the  $^{\mathrm{K}}$ DMDGP have power of two solutions. Yet, the computational trend remained unexplained. The main contribution of this paper is a proof that the set of YES instances of the  $^{\mathrm{K}}$ DMDGP such that  $|X|$  is a power of two has Lebesgue measure 1 in the set of all YES instances of the  $^{\mathrm{K}}$ DMDGP. The statement is based on the assumption that we consider edge weight functions whose range consists of real numbers. We also partially extend this result to graphs which are more realistic protein models.

# 3.1. Statements holding with probability 1

In the following, we assume that the probability of any point of  $\mathbb{R}^K$  belonging to any given subset of  $\mathbb{R}^K$  having Lebesgue measure zero is equal to zero. Based on this assumption, when we state “ $(\forall p \in P F(p))$  with probability 1” for a certain well-formed formula  $F$  with a free variable ranging over a set  $P$  having a strictly positive Lebesgue measure, we really mean that there exists a Lebesgue measurable subset  $Q \subseteq P$ , with Lebesgue measure equal to that of  $P$ , such that  $\forall p \in Q F(p)$ . Equivalently, statements holding with probability 1 should be taken to mean that the set of  $^K$ DMDGP instances and partial realizations  $x$  for which the statements do not hold has Lebesgue measure zero in the set of all  $^K$ DMDGP instances and partial realizations.

In this paper, zero Lebesgue measure sets are associated to cases that occur whenever pairs of real numbers (such as components of vectors) happen to be equal: we remark that affine subspaces of  $\mathbb{R}^K$  defined by linear equations all have zero Lebesgue measure in  $\mathbb{R}^K$ . For example, the set of all pairs of points  $y, y' \in \mathbb{R}^K$  having  $t$  equal components (say, the first  $t$ ), is described by the  $t$  linear equations  $\forall j \leq t$  ( $y_j = y_j'$ ) and therefore has zero Lebesgue measure in the set of all pairs of points in  $\mathbb{R}^K$ . In other words, when uniformly sampling pairs of points randomly from  $\mathbb{R}^K$ , the event that they should end up having at least one equal component has probability zero. This is similar to the spirit of Graver's definition of genericity using the nontrivial minors of the complete rigidity matrix [13]. The guiding principle in most of the results in this paper is that whenever a logical case of a proof requires two (or more) real values, sampled in a uniform distribution over a set of positive Lebesgue measure, to be equal, then this case has probability 0 to occur. The sampling usually refers to choosing instances randomly from the  $^K$ DMDGP set, and more precisely to the real weights assigned to the edges of  $G$ . This directly translates to a uniform sampling over vectors of  $\mathbb{R}^K$  (the components of the realizations associated to each instance).

We remark that the notion we consider is different from the usual genericity notion employed in rigidity theory [6], which requires distances to be algebraically independent over  $\mathbb{Q}$ . In order to solve the  $^{\mathrm{K}}$  DMDGP practically, we deal with instances whose edge weights are encoded in floating points, which are certainly not algebraically independent over  $\mathbb{Q}$  (since they are themselves a subset of  $\mathbb{Q}$ ). For example, Lemma 4.2 would not hold under algebraic independence (the intersection of  $K + 1$  "generic spheres" in  $\mathbb{R}^K$  is empty), but it holds under our weaker requirement.

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

# 4. Formal definitions: the $^{K}$DMDGP and the BP algorithm

For a set $U = \{x_{i} \in \mathbb{R}^{K} \mid i \leq K + 1\}$ of points in $\mathbb{R}^{K}$, let $D$ be the symmetric matrix whose $(i, j)$-th component is $\| x_{i} - x_{j} \|^{2}$ for all $i, j \leq K + 1$ and let $D'$ be $D$ bordered by a left $(0, 1, \ldots, 1)^{\top}$ column and a top $(0, 1, \ldots, 1)$ row (both of size $K + 2$). Then the Cayley-Menger formula states that the volume $\Delta_{K}(U)$ of the $K$-simplex on $U$ is given by $\Delta_{K}(U) = \sqrt{\frac{(-1)^{K + 1}}{2^{K}(K')^{2}}} |D'|$. The strict simplex inequalities are given by $\Delta_{K}(U) &gt; 0$. For $K = 2$, these reduce to strict triangle inequalities. We remark that only the distances of the simplex edges are necessary to compute $\Delta_{K}(U)$, rather than the actual points in $U$; the needed information can be encoded as a $(K + 1)$-clique with these distances as edge weights.

Let $n = |V|$ and $m = |E|$. For all $v \in V$, let $N(v) = \{u \in V \mid \{u, v\} \in E\}$ be the star of vertices around $v$ (also called the adjacencies of $v$); for a directed graph $(V, A)$, where $A \subseteq V \times V$, we denote the outgoing star by $N^{+}(v) = \{u \in V \mid (v, u) \in A\}$. For an order $&lt;$ on $V$, let $\gamma(v) = \{u \in V \mid u &lt; v\}$ be the set of predecessors of $v$, and let $\rho(v) = |\gamma(v)| + 1$ be the rank of $v$ in $&lt;$. For $V' \subseteq V$, we denote by $G[V']$ the subgraph of $G$ induced by $V'$. For a finite set $M$, let $\mathcal{P}(M)$ be its power set. For a sequence $x = (x_{1}, \ldots, x_{n})$ and a subset $U \subseteq \{1, \ldots, n\}$ we let $x[U]$ be the subsequence of $x$ indexed by $U$. If $x$ is an initial subsequence of $y$, then $y$ is an extension of $x$. We denote a clique on $q$ vertices by $\mathbf{K}_q$.

GENERALIZED DISCRETIZABLE MOLECULAR DISTANCE GEOMETRY PROBLEM ($^{K}$DMDGP). Given an undirected graph $G = (V, E)$, an edge weight function $d: E \to \mathbb{R}_{+}$, an integer $K &gt; 0$, a subset $V_{0} \subseteq V$ with $|V_{0}| = K$, a partial realization $\widehat{x}: V_{0} \to \mathbb{R}^{K}$ valid for $G[V_{0}]$, and a total order $&lt;$ on $V$ such that for each $v$ with $\rho(v) &gt; K$ there is a set $U_{v} \subseteq N(v) \cap \gamma(v)$ with the following properties:

$$
\{v \in V \mid \rho (v) \leq K \} = V _ {0}; \tag{2}
$$

$$
\forall u \in U _ {v} \quad (\rho (v) - K \leq \rho (u) \leq \rho (v) - 1); \tag{3}
$$

$$
\forall v \in V \setminus V _ {0} \quad (G [ U _ {v} ] = \mathbf {K} _ {K} \wedge \Delta_ {K - 1} (U _ {v}) &gt; 0), \tag{4}
$$

decide whether there is a valid extension $x:V\to \mathbb{R}^K$ of $\widehat{x}$

Condition (2) requires the first $K$ vertices to induce a $K$-clique in $G$. Condition (3) requires the $K$ immediate predecessors of $v$ to be also adjacent to it. Condition (4), which was not mentioned in the informal discussion in Sections 1 and 2, prevents the realizations $x$ of $G$ from being flexible because of an algebraic dependence on the components of $x$. Specifically, requiring the Cayley-Menger determinants to be strictly positive implies that no $K$-clique determined by immediate predecessor will be embedded as a volume 0 simplex in $\mathbb{R}^{K-1}$ (Fig. 3 shows the effect of a triangle embedded in a line instead of in a plane: the SIP fails to hold).

Two further remarks are in order.

- Whenever vertex $v$ is being considered in the BP algorithm, all its predecessors have already been placed. Hence, all of the distances between all predecessors are already known; thus, the BP can also solve instances for which $G[U_v]$ is not the full $K$-clique, although they are not formally in the $^K$DMDGP.
- Edge weights are real numbers in $^K$DMDGP instances: this naturally makes the Turing machine model difficult to apply. We observe, however, that it is not known whether the DGP is in NP, since the certifying realizations might have irrational components even in the case when the weights are restricted to be integer (just take an equilateral triangle with unit weights). In practice, our algorithms work with floating point numbers, so whenever we say "exact" we mean "exact in theory". In practice, solutions will be floating point approximations of the exact solutions. This, however, holds for all existing methods targeting DGP-type problems.

## 4.1. Branching and pruning

Let $G$ be a $^K$DMDGP instance. Consider $v \in V$ with rank $\rho(v) = i &gt; K$, let $G^v = G[\gamma(v) \cup \{v\}]$ and $x$ be a valid realization of $G[\gamma(v)]$. We characterize the number of extensions of $x$ valid for $G^v$ in the following lemmata (which also hold for the DDGP). Lemma 4.1 (resp. 4.2) essentially state that, under the given conditions, $G[\{v\} \cup (N(v) \cap \gamma(v))]$ is a rigid (resp. uniquely rigid) graph. The results in this section are not new, but we list the proofs here because these lemmata form the basic toolbox for what is to follow.

**Lemma 4.1.** If $|N(v) \cap \gamma(v)| = K$ then there are at most two distinct extensions of $x$ that are valid for $G^v$. If one valid extension exists, then with probability 1 there are exactly two distinct valid extensions.

**Proof.** Since $|N(v) \cap \gamma(v)| = K$, $U_v = N(v) \cap \gamma(v)$ and $v$ is at the intersection of exactly $K$ spheres in $\mathbb{R}^K$ (each centered at $x_u$ with radius $d_{uv}$, where $u \in U_v$). The position $z \in \mathbb{R}^K$ of $v$ must then satisfy:

$$
\forall u \in U _ {v} \quad \| z - x _ {u} \| = d _ {u v} \Rightarrow \| z \| ^ {2} - 2 x _ {u} \cdot z + \| x _ {u} \| ^ {2} = d _ {u v} ^ {2}. \tag{5}
$$

As in [11], we choose an arbitrary $w \in U_v$, say $w = \max_{&lt;} U_v$, and subtract from the row of Eq. (5) indexed by $w$ the other rows of (5), obtaining the system:

$$
\left. \begin{array}{l}
\forall u \in U _ {v} \setminus \{w \} \quad 2 \left(x _ {u} - x _ {w}\right) \cdot z = \left(\| x _ {u} \| ^ {2} - d _ {u v} ^ {2}\right) - \left(\| x _ {w} \| ^ {2} - d _ {w v} ^ {2}\right) \\
\| z \| ^ {2} - 2 x _ {w} \cdot z + \| x _ {w} \| ^ {2} = d _ {w v} ^ {2}
\end{array} \right\}. \tag{6}
$$

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

System (6) consists of a set of $K - 1$ linear equations and a single quadratic equation in the $K$-vector $z$. We write the linear equations as the system $Az = b$, where the $(u, j)$-th component of $A$ is $2(x_{uj} - x_{wj})$, the $u$-th component of $b$ is $\|x_u\|^2 - \|x_w\|^2 - d_{uv}^2 + d_{uvv}^2$. $A$ is a $(K - 1) \times K$ matrix and $b \in \mathbb{R}^{K-1}$. By the strict simplex inequalities, $A$ has full rank (for otherwise the linear dependence condition $\sum_{u \neq w} \xi_u(x_u - x_w) = 0$, for some coefficients $\xi_u$, implies that $x_w$ is in the span of $\{x_u \mid u \in U_v\}$, and hence that $\Delta_{K-1}(U_v) = 0$); so without loss of generality assume that the square matrix $B$ formed by the first $K-1$ columns of $A$ is invertible. Let $z_B$ be the vector consisting of the first $K-1$ components of $z$; then the linear part (first $K-1$ equations) of (6) yields $z_B = B^{-1}(b - Nz_K)$ as a function $z_B(z_K)$ of $z_K$, where $N = 2(x_{uK} - x_{wK} \mid u \in U_v \setminus \{w\}) \in \mathbb{R}^{K-1}$. After replacement of $z_B$ in (6) with $z_B(z_K)$, we obtain the following quadratic equation in $z_K$:

$$
\left(\| \tilde {N} \| ^ {2} + 1\right) z _ {K} ^ {2} - 2 \left(\left(\tilde {b} + x _ {w B}\right) \tilde {N} + x _ {w K}\right) z _ {k} + \left(\| x _ {w B} - \tilde {b} \| ^ {2} + x _ {w K} ^ {2} - d _ {u v} ^ {2}\right) = 0, \tag {7}
$$

where $\tilde{b} = B^{-1}b$ and $\tilde{N} = B^{-1}N$. If the discriminant of (7) is negative then no extension of $\tilde{x}$ to $v$ is possible and the result follows. If the discriminant is nonnegative, (7) has solutions $z_{K}^{\prime}, z_{K}^{\prime \prime}$ yielding points $z^{\prime} = (z_{B}(z_{K}^{\prime}), z_{K}^{\prime})$ and $z^{\prime \prime} = (z_{B}(z_{K}^{\prime \prime}), z_{K}^{\prime \prime}) \in \mathbb{R}^{K}$, which are distinct with probability 1 because the discriminant is zero with probability 0. The extended realizations, distinct with probability 1, are given by $(x, z^{\prime})$ and $(x, z^{\prime \prime})$.

Lemma 4.2. If $|N(v) \cap \gamma(v)| &gt; K$ then, with probability 1, there is at most one extension of $x$.

Proof. Consider a subset $S \subseteq N(v) \cap \gamma(v)$ such that $|S| = K + 1$ and $S \supseteq U_v$. Either there is at least one point $x_v$ such that $(x, x_v)$ is a realization of $G[S \cup \{v\}]$ that is valid w.r.t. the system:

$$
\forall u \in S \quad \sum_ {k \leq K} \left(x _ {v k} ^ {2} - 2 x _ {u k} x _ {v k} + x _ {u k} ^ {2}\right) = d _ {u v} ^ {2}, \tag {8}
$$

or the system has no solution. In the latter case, the result follows, so we assume now that there is a point $x_{v}$ satisfying (8). Since the points $x_{u}$ are known for all $u \in S$, (8) is a quadratic system with $K$ variables and $K + 1$ equations. As in the proof of Lemma 4.1, we derive an equivalent linear system from (8). Since $d$ satisfies the strict simplex inequalities on $U_{v}$ with probability 1 and $S \supseteq U_{v}$, by [7] $\{x_{u} \mid u \in S\}$ are not co-planar and the system has exactly one solution, as claimed.

Lemma 4.3. With the notation of Lemma 4.1, if $\tilde{x}$ is a valid realization for $G[U_v]$, then $z''$ is a reflection of $z'$ with respect to the hyperplane through the $K$ points of $\tilde{x}$.

Proof. Any sphere in $\mathbb{R}^K$ is symmetric with respect to any hyperplane through its center; so the intersection of up to $K$ spheres in $\mathbb{R}^K$ is symmetric with respect to the hyperplane containing all the centers.

Remark 4.4. Reflections with respect to hyperplanes are isometries, and can therefore be represented by linear operators. If $a \in \mathbb{R}^K$ is the unit normal vector to a hyperplane $H$ containing the origin, then the reflection operator $R_0$ w.r.t. $H$ can be expressed in function of the standard basis by the matrix $I - 2aa^\top$, where $I$ is the $K \times K$ identity matrix [5]. Let $H$ be a hyperplane with equation $a^\top x = a_0$ (with $a_0 \neq 0$) and $a_i$, for some $1 \leq i \leq K$, be the nonzero coefficient of smallest index in $a$. Then, the reflection operator $R$ acting on a point $p \in \mathbb{R}^K$ w.r.t. $H$ is given by $R(p) = R_0(p - \frac{a_0}{a_i} e_i) + \frac{a_0}{a_i} e_i$, where $e_i \in \mathbb{R}^K$ is the unit vector with 1 at index $i$ and 0 elsewhere: we first we translate $p$ so that we can reflect it using $R_0$ w.r.t. the translation of $H$ containing the origin, then we perform the inverse translation of the reflection.

## 4.2. The BP search tree

We denote the BP binary search tree by $\mathcal{T} = (\mathcal{V},\mathcal{A})$. The tree is directed from the root to the leaf nodes, which are triplets $\alpha = (x(\alpha),\lambda (\alpha),\mu (\alpha))$. Informally, $x(\alpha)$ is a realization from the root to node $\alpha$, $\lambda (\alpha)$ is 0 or 1 according as to whether $\alpha$ is a left or right subnode of its parent, and $\mu$ is $\boxplus$ if $x(\alpha)$ is feasible and $\boxminus$ otherwise.

More formally, for $\alpha \in \mathcal{T}$ we denote by $\rho(\alpha)$ the unique path from the root node $r$ of $\mathcal{T}$ to $\alpha$, and by $\alpha^{-}$ the unique parent node of $\alpha$ (unless $\alpha = r$, in which case we define $r^{-} = r$). The symbol $x(\alpha)$ is defined recursively to denote an extension of the realization $x^{-} = x(\alpha^{-})$ found on $\rho(\alpha^{-})$. The symbol $\lambda(\alpha) \in \{0, 1\}$ distinguishes whether $\alpha$ is a "left" or a "right" subnode of $\alpha^{-}$. More precisely, let $\alpha$ be a node at level $i$ in $\mathcal{T}$, $v = \rho^{-1}(i)$, $\tilde{x}$ be a partial realization of $G[U_v]$, and $a_v^\top x = a_{v0}$ be the equation of the hyperplane through the points of $\tilde{x}$, which is $(K - 1)$-dimensional by (4). Assuming that $u = \rho^{-1}(i - 1)$ and $a_v \in \mathbb{R}^K$ is oriented so that $a_v \cdot a_u \geq 0$; then:

$$
\lambda (\alpha) = \left\{ \begin{array}{l l} 0 &amp; \text{if } a _ {v} ^ {\top} x (\alpha) _ {i} \leq a _ {v 0} \\ 1 &amp; \text{if } a _ {v} ^ {\top} x (\alpha) _ {i} &gt; a _ {v 0}. \end{array} \right. \tag {9}
$$

Lastly, $\mu(\alpha) = \boxplus$ if $x(\alpha)$ is a valid extension of $x^{-}$, in which case the node is said to be feasible, and $\mu = \boxminus$ otherwise. This allows us to retrieve the set $X$ of all valid realizations of $G$ by simply traversing $\mathcal{T}$ backwards from the leaf nodes marked $\boxplus$ up to $r$.

Although the notation we introduced to describe the BP algorithm looks overly formal and complicated, it allows us to give rigorous proofs of subsequent material.

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

Algorithm 1 The Branch and Prune algorithm.
Require: Partial realization $\tilde{x}$ of first $K$ vertices of $G$
Ensure: Set $X$ of valid realizations of $G$
1: Let $\alpha = (\tilde{x}_1,0,\boxplus)$ and $\alpha^{\prime} = (\tilde{x}_{1},1,\boxminus)$
2: Initialize $\mathcal{V} = \{\alpha ,\alpha^{\prime}\}$ and $\mathcal{A} = \{(r,\alpha),(r,\alpha^{\prime})\}$
3: for $1 &lt; i\leq K$ do
4: Let $\alpha = (\tilde{x}_i,0,\boxplus),\alpha^{\prime} = (\tilde{x}_{i},1,\boxminus),\beta = (\tilde{x}_{i - 1},0,\boxplus)$
5: Let $\mathcal{V}\gets \mathcal{V}\cup \{\alpha ,\alpha '\}$ and $\mathcal{A}\gets \mathcal{A}\cup \{(\beta ,\alpha),(\beta ,\alpha^{\prime})\}$
6: end for
7: BRANCHANDPRUNE $(K + 1,(\tilde{x}_K,0,\boxplus))$
8: Let $X = \{x(\theta)\mid \theta \in \mathcal{V}\land |N^{+}(\theta)| = 0\land \mu (\theta) = \boxplus \}$
9: stop
10:
function BRANCHANDPRUNE(i,β):
if $i &gt; n\lor \mu = \boxminus$ then
11: return
12: end if
13: Let $v = \rho^{-1}(i)$
14: Compute the equation $a_v^\top x = a_{v0}$ of the hyperplane through $x[U_v]$
15: Let $Z = \{z^{\prime},z^{\prime \prime}\}$ be extensions of $x(\beta)$ to $v$, and $Z^{\prime} = Z$
16: for $z\in Z$ do
17: if $\exists \{u,v\} \in E\| x(\beta)_u - z\| \neq d_{uv}$ then
18: Let $Z = Z\setminus \{z\}$
19: end if
20: end for
21: if $Z = \{z^{\prime},z^{\prime \prime}\}$ then
22: if $a_v^\top z' \leq a_{v0}$ then
23: Let $\alpha = (z',0,\boxplus),\alpha' = (z'',1,\boxplus)$
24: else
25: Let $\alpha = (z'',0,\boxplus),\alpha' = (z',1,\boxplus)$
26: end if
27: else if $Z = \{z\}$ then
28: if $a_v^\top z \leq a_{v0}$ then
29: Let $\alpha = (z,0,\boxplus),\alpha' = (Z' \setminus \{z\},1,\boxplus)$
30: else
31: Let $\alpha = (z,1,\boxplus),\alpha' = (Z' \setminus \{z\},0,\boxplus)$
32: end if
33: else
34: return
35: end if
36: Let $\mathcal{V}\gets \mathcal{V}\cup \{\alpha ,\alpha '\}$ and $\mathcal{A}\gets \mathcal{A}\cup \{(\beta ,\alpha),(\beta ,\alpha^{\prime})\}$
37: for $\theta \in N^{+}(\beta)$ such that $\mu (\theta) = \boxplus$ do
38: BRANCHANDPRUNE(i+1,θ)
39: end for
40: return

We remark that Algorithm 1 differs from the original BP formulation [28] because it applies to $K$ dimensions and explicitly stores several details of the binary search tree.

Lemma 4.5. At termination of Algorithm 1, $X$ contains all valid realizations of $G$ extending $\tilde{x}$.

Proof. $Z$ (in Step 17) exists with probability 1 by Lemma 4.1. Every realization in $X$ is valid because of Steps 17 and 19-20. No other valid extension of $\tilde{x}$ exists because of Lemmata 4.1-4.2. $\square$

The realizations in $X$ are incongruent apart perhaps from between 1 and $n - K$ reflections along the hyperplane defined by $\tilde{x}_1, \ldots, \tilde{x}_v$. The exact number of these reflections depends on the order on $V$ and the edges in $G$.

We now partition $\mathcal{V}$ in pairwise disjoint subsets $\mathcal{V}_1, \ldots, \mathcal{V}_n$ where for all $i \leq n$ the set $\mathcal{V}_i$ contains all the nodes of $\mathcal{V}$ at level $i$ of the tree $\mathcal{T}$. We show in Proposition 4.6 that no level of $\mathcal{T}$ has two distinct feasible nodes having respectively one and two feasible subnodes.

Proposition 4.6. With probability 1, there is no level $i \leq n$ having two distinct feasible nodes $\beta, \theta \in \mathcal{V}_i$ such that $|\{\alpha \in N^+(\beta) \mid \mu(\alpha) = \boxplus\}| = 1$ and $|\{\alpha \in N^+(\theta) \mid \mu(\alpha) = \boxplus\}| = 2$.

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

Proof. We show that for all $i \leq n$ the event of having two distinct nodes $\beta, \theta \in \mathcal{V}_i$, with $\rho^{-1}(i) = v$, such that $\beta$ has one feasible subnode and $\theta$ has two has probability 0. Consider $T_v = N(v) \cap \gamma(v)$: if $|T_v| = K$ then, by Lemma 4.1, $\beta$ should have exactly two feasible subnodes with probability 1. Since by hypothesis it only has one, the event $|T_v| = K$ occurs with probability 0. Since $|T_v| \geq K$ by (4), the event $|T_v| &gt; K$ occurs with probability 1. Thus by Lemma 4.2 there is, with probability 1, at most one valid realization extending the partial realization at $v$, which means that the two feasible subnodes of $\theta$ represent the same realization, an event that occurs with probability 0.

We remark that Proposition 4.6 also holds for the DDGP, provided $U_v$ is chosen in Algorithm 1 as any subset of $N(v) \cap \gamma(v)$ satisfying the constraints of Eq. (4).

Lastly, we emphasize the fact that for all $\ell \in \{i,\dots ,n\}$ and for all $\alpha \in \mathcal{V}_{\ell}$ the set $\mathsf{p}(\alpha)\cap \mathcal{V}_i$ has a unique element, as it contains the unique node at level $i$ on the path from $\alpha$ to the BP tree root node.

# 5. Geometry in BP trees

The most important result of this section is that, for any valid realization $y \in X$, if the BP tree branches at level $i = \rho(v)$ on the path to $y$ and both branches continue to the last level, then the realization obtained by reflecting all the points of $y$ past the $(i-1)$-st vertex through the hyperplane defined by $y[U_v]$ is also valid with probability 1. We remark that the results in this section only apply to the $^{\kappa}$ DMDGP (not to the DDGP, as shown in the counterexample of Fig. 8).

We need to emphasize those BP branchings which carry on to feasible leaf nodes along both branches. For $y \in X$ and a vertex $v \in V \setminus V_0$ we denote $\Upsilon(y, v)$ the following property:

$\Upsilon(y, v)$: there are feasible leaf nodes $\beta, \beta' \in \mathcal{V}_n$ such that $x(\beta) = y$, $\mathsf{p}(\beta) \cap \mathcal{V}_{\rho(v)-1} = \mathsf{p}(\beta') \cap \mathcal{V}_{\rho(v)-1}$ and $\mathsf{p}(\beta) \cap \mathcal{V}_{\rho(v)} \neq \mathsf{p}(\beta') \cap \mathcal{V}_{\rho(v)}$.

In other words, $\Upsilon(y, v)$ holds whenever $\mathsf{p}(\beta) \cap \mathcal{V}_{\rho(v)-1}$ contains a feasible node with two feasible subnodes leading to different valid realizations.

# 5.1. Partial reflection operators

With $\Upsilon(y, v)$ true, we let $R^v$ be the Euclidean reflection operator with respect to the hyperplane through $y[U_v]$ (as discussed in Remark 4.4). Define

$$
\tilde {R} ^ {v} = I ^ {\rho (v) - 1} \times \left(R ^ {v}\right) ^ {n - \rho (v)},
$$

i.e., $\tilde{R}^v y = (y_1, \ldots, y_{\rho(v) - 1}, R^v y_{\rho(v)}, \ldots, R^v y_n)$ for any realization $y$. This is a partial reflection of $y$ which only acts on vertices past rank $\rho(v) - 1$.

We remark that the matrix representing $R^v$ could change depending on $y$. Since we wish $R^v$ to represent the effect of a reflection at level $v$ rather than the reflection matrix itself, we introduce the following technicality. Consider an equivalence relation on the set of all possible (not necessarily feasible) realizations $V \to \mathbb{R}^K$ extending $\tilde{x}$. Let $E_{\mathrm{d}} \subseteq E$ be the set of edges $\{u, v\} \in E$ such that $|\rho(v) - \rho(u)| \leq K$, and $G_{\mathrm{d}}$ be the subgraph of $G$ defined by $E_{\mathrm{d}}$. Two realizations $y, y'$ are equivalent if both are valid for $G_{\mathrm{d}}$. Now $R^v$ can be formally defined as the set of matrices representing $R^v$ over all realizations of the same equivalence class. This definition carries over to $\tilde{R}^v$.

The following is an easy corollary to Lemma 4.3, and states, with the formal notation introduced for Algorithm 1, that two feasible subnodes of a feasible node are associated with partial realizations whose last components are $K$-vectors which are reflections of each other.

Corollary 5.1. Let $\alpha \in \mathcal{V}_{i-1}$ for some $i &gt; 1$, $v = \rho^{-1}(i)$ and $N^{+}(\alpha) = \{\eta, \beta\}$ with $\mu(\eta) = \mu(\beta) = \boxplus$. Then $x(\eta)_v = R^v x(\beta)_v$.

# 5.2. Probability of conditional events

In most subsequent results, we assume the considered $^{\kappa}$ DMDGP instance to be a YES one, so that probabilities are conditional to this event.

The following remark is in order.

Remark 5.2. If $\Upsilon(y, v)$ holds for some $y \in X$ and $v \in V \setminus V_0$, then by definition there are feasible leaf nodes in the BP tree, which implies that the considered $^{\kappa}$ DMDGP instance is YES.

An important consequence of Remark 5.2 is that all statements assuming $\Upsilon(y, v)$ and claiming a result with probability 1 implicitly also assume that the probability is conditional to the event of the $^{\kappa}$ DMDGP instance being a YES one. In particular, since the instance is YES, certain points must be placed at certain distances with probability 1 (notably, at distances satisfying the Eqs. (1)), for otherwise the instance would be NO. This is evident in Proposition 5.4, Corollaries 5.6, 5.7, and Theorem 5.9, where we state that certain real scalars and vectors must belong to certain finite sets with probability 1. The sense of these assertions, in this context, is that the Lebesgue measure of the set of YES instances not satisfying the result is zero in the set of all YES instances.

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

![img-3.jpeg](img-3.jpeg)
(a)  $a_{v}^{0\top}x(\eta)_{v} &gt; a_{v0}^{0}$  and  $a_{v}^{1\top}x(\beta)_{v} &gt; a_{v0}^{1}$ .

![img-4.jpeg](img-4.jpeg)
(b)  $a_{v}^{0\top}x(\eta)_{v} &gt; a_{v0}^{0}$  and  $a_{v}^{1\top}x(\beta)_{v} &lt; a_{v0}^{1}$ .
Fig. 4. Proof of Lemma 5.3: Case 4(a) shows the contradiction deriving from  $\lambda (\eta) = \lambda (\beta) = 0$  (or  $x(\beta)_v\neq R^u x(\eta)_v$ ), and case 4(b) the situation that actually occurs.

# 5.3. Reflections in the BP tree

We build towards the main result of the section; we start with a technical lemma which relates the position of reflected points in the realizations with the "left/right"  $\lambda$ -components of the corresponding BP node triplets.

Lemma 5.3. Let  $\alpha \in \mathcal{V}_{i-1}$  for some  $i &gt; 1$  such that  $N^{+}(\alpha) = \{\eta', \beta'\}$ , let  $u = \rho^{-1}(i)$ , let  $v$  with  $\rho(v) = \ell$  be such that  $\ell &gt; i$ , and consider two feasible nodes  $\eta, \beta \in \mathcal{V}_{\ell}$  such that  $\{\eta'\} = \mathsf{p}(\eta) \cap \mathcal{V}_i$  and  $\{\beta'\} = \mathsf{p}(\beta) \cap \mathcal{V}_i$ . Further,  $\forall i \leq j \leq \ell$  let  $w = \rho^{-1}(j)$ ,  $\mathsf{p}(\eta) \cap \mathcal{V}_j = \{\eta''\}$  and  $\mathsf{p}(\beta) \cap \mathcal{V}_j = \{\beta''\}$ . Then, with probability 1, the following statements are equivalent:

(i)  $\forall i\leq j\leq \ell x(\beta^{\prime \prime})_w = R^u x(\eta^{\prime \prime})_w;$
(ii)  $\forall i\leq j\leq \ell \lambda (\eta^{\prime \prime}) = 1 - \lambda (\beta^{\prime \prime}).$

Proof. Let  $a_v^{0\top}x = a_{v0}^0, a_v^{1\top}x = a_{v0}^1$  be the equations of the hyperplanes  $H_\eta, H_\beta$  defined respectively by  $x(\eta)[U_v]$  and  $x(\beta)[U_v]$ , with the normals oriented as explained in Remark 4.4. We prove by induction on  $\ell - i$  that the following statement is equivalent to (i) and (ii):

(iii) for all  $i \leq j \leq \ell$ ,  $x(\beta'')_w = R^u x(\eta'')_w$  and  $a_u \cdot a_w^0 = a_u \cdot a_w^1$ , where  $a_w^0$  and  $a_w^1$  are the normal vectors of the hyperplanes  $H_{\eta''}$  and  $H_{\beta''}$  oriented as usual.

The induction starts because, if  $\ell = i$ , then (i)-(iii) hold simultaneously. Indeed,  $\eta = \eta'$  and  $\beta = \beta'$ , hence  $x(\beta)_v = R^u x(\eta)_v$  (Lemma 4.3) and  $\lambda(\eta) = 1 - \lambda(\beta)$  (Algorithm 1, Steps 25 and 27). In addition, we have  $H_\eta = R^u H_\beta$ , therefore  $|a_u \cdot a_v^0| = |a_u \cdot a_v^1|$ . Because the orientation of  $a_v^0$ ,  $a_v^1$  is such that  $a_u \cdot a_v^0$ ,  $a_u \cdot a_v^1 \geq 0$ , the result holds. Assume that the equivalence stated above holds for level  $\ell - 1$ , we show that it is still the case at level  $\ell$ . In the sequel, denote  $t = \rho^{-1}(\ell - 1)$ .

$(\mathrm{iii})\Rightarrow (\mathrm{i})$  : trivial.

$(\mathrm{i}) \Rightarrow (\mathrm{ii})$ . Suppose for all  $i \leq j &lt; \ell$ ,  $x(\beta'')_w = R^u x(\eta'')_w$  and  $\lambda(\eta'') = 1 - \lambda(\beta'')$  (by the induction hypothesis, both statements are equivalent). Hence,  $H_{\eta''} = R^u H_{\beta''}$  holds for all  $j$ , because the  $K$  points generating the hyperplanes either belong to  $H_\alpha$ , or are reflections of each other. This is true in particular if we choose  $\eta'', \beta'' \in \mathcal{V}_{\ell-1}$ . In addition, if we use the induction hypothesis  $(\mathrm{i}) \Rightarrow (\mathrm{iii})$ , we have  $a_u \cdot a_\ell^0 = a_u \cdot a_\ell^1$ , so  $a_\ell^0, a_\ell^1$  are directed similarly w.r.t.  $a_u$ , and hence  $\lambda(\eta) = 1 - \lambda(\beta)$  if and only if  $x(\beta)_v = R^u x(\eta)_v$  (see Fig. 4).
(ii)  $\Rightarrow$  (iii). Suppose for all  $i \leq j \leq \ell$ ,  $\lambda(\eta'') = 1 - \lambda(\beta'')$ . By the previous result, we also know that for all  $i \leq j \leq \ell$ ,  $x(\beta'')_w = R^u x(\eta'')_w$ . It remains to prove that  $a_u \cdot a_v^0 = a_u \cdot a_v^1$ , i.e. that the angles  $\theta_v^0$  and  $\theta_v^1$  formed by these vectors have the same cosine. Notice once again that  $H_\eta = R^u H_\beta$ . By induction, we know that the angles  $\theta_\ell^0$ ,  $\theta_\ell^1$  formed by  $a_u$  and respectively  $a_\ell^0$ ,  $a_\ell^1$ , have the same cosine. With probability 1, the hyperplanes  $H_\eta$ ,  $H_\beta$  are not parallel, hence their normal vectors cannot be identical, therefore,  $\theta_\ell^0 = -\theta_\ell^1$  (see the illustration on Fig. 5). Denote  $\theta^0$ ,  $\theta^1$  the angles formed respectively by  $a_\ell^0$  and  $a_v^0$ , and by  $a_\ell^1$  and  $a_v^1$ . We also have,  $H_{\eta''} = R^u H_{\beta''}$ , where  $\eta'', \beta'' \in \mathcal{V}_{\ell-1}$ , hence the normal vectors of these 4 hyperplanes are also symmetric, which implies  $\theta^0 = -\theta^1$  or  $\theta^0 = \pi - \theta^1$ . By the definition of  $a_v^0$  and  $a_v^1$  (page remark), since the scalar products are positive,  $-\pi/2 \leq \theta^0$ ,  $\theta^1 \leq \pi/2$ , thus  $\theta^0 = -\theta^1$ . Therefore,  $\theta_v^0 = \theta_\ell^0 + \theta^0 = -\theta_\ell^1 - \theta^1 = -\theta_v^1$ , which concludes this part of the proof.

The following important proposition considers vertices  $u, v \in V$  such that  $|\rho(v) - \rho(u)| = K + 1$  and states that there are only two possible values of  $\| y_u - y_v \|$  if  $y$  is to be a feasible realization.

Proposition 5.4. Consider a subtree  $\mathcal{T}'$  of  $\mathcal{T}$  consisting of  $K + 2$  consecutive levels  $i - K - 1, \ldots, i$  (where  $i \geq 2K + 1$ ), rooted at a single node  $\eta$  and such that all nodes at all levels are marked. Let  $p = 2^{K + 1}$  and consider the set  $Y' = \{y^j \mid j \leq p\}$  of partial

---

L. Liberti et al./ Discrete Applied Mathematics 165 (2014) 213-232

![img-5.jpeg](img-5.jpeg)
Fig. 5. Proof of Lemma 5.3: illustration of the fact that  $a_{u} \cdot a_{v}^{0} = a_{u} \cdot a_{v}^{1}$ .

![img-6.jpeg](img-6.jpeg)
Fig. 6. Proof of Proposition 5.4 in  $\mathbb{R}^2$ . The arrangement of three segments gives rise, in general, to two distances  $r$ ,  $r'$  between root and leaves.

realizations of  $G$  at the leaf nodes  $\{\alpha_j \mid j \leq p\}$  of  $\mathcal{T}'$ , such that  $\forall j \leq p (y^j = x(\alpha_j))$ . Let  $u = \rho^{-1}(i - K - 1)$  and  $v = \rho^{-1}(i)$ . Then with probability 1 there are two distinct positive reals  $r, r'$  such that  $\| y_u^j - y_v^j \| \in \{r, r'\} = H^{uv}$  for all  $j \leq p$ .

Proof. Fig. 6 shows a graphical proof sketch for  $K = 2$ . With a slight abuse of notation, for a vertex  $w \in V$  in this proof we denote by  $R^w$  the set of all reflections at level  $w$ . We order the  $\alpha_j$  nodes (and the corresponding  $y^j$ ) so that the action of  $R^v$  on  $(\alpha_1, \ldots, \alpha_p)$  is the permutation  $\prod_{j \bmod 2=1} (j, j+1)$ . Let  $t = \rho^{-1}(i-1)$ . Since all nodes are feasible,  $\| y_v^j - y_t^j \| = d_{tv}$  and  $\| y_u^j - y_t^j \| = d_{ot}$  for all  $j \leq p$  (we remark that  $\{t, v\}$  and  $\{u, t\}$  must be in  $E$  by the  $^K$ DMDGP definition). With probability 1, the segments through  $y_u^j$  and  $y_t^j$  (where  $j \leq p$ ) do not respectively lie within the hyperplanes defining the reflections  $R^v$ ; and the same holds for the segments through  $y_t^j$  and  $y_v^j$ . Thus, there is a set  $Q$  of positive reals  $r_1, \ldots, r_p$  s.t. for all  $j \leq p$  with  $j \bmod 2 = 1$  we have  $\| y_u^j - y_v^j \| = r_j$  and  $\| y_u^{j+1} - y_v^{j+1} \| = r_{j+1}$ , which shows  $|Q| \leq p = 2^{K+1}$ . By Lemma 5.3, the action of  $R^i$  on  $(\alpha_1, \ldots, \alpha_p)$  is the permutation  $\prod_{j \bmod 4=1} (j, j+3)(j+1, j+2)$ : this implies that  $r_j = r_{j+3}$  and  $r_{j+1} = r_{j+2}$  for all  $j \bmod 4 = 1$ , which shows  $|Q| \leq p/2 = 2^K$ . Inductively, for a vertex  $w$  s.t.  $i - K \leq \rho(w) \leq i - 1$  the action of  $R^w$  is  $\prod_{j \bmod 2^{j-\rho(w)+1}} (j, j+2^{i-\rho(w)+1}-1)(j+1, j+2^{i-\rho(w)+1}-2) \cdots (j+2^{i-\rho(w)}-1, j+2^{i-\rho(w)})$ , which implies that  $|Q| \leq 2^{K+1-i+\rho(w)}$ . Therefore  $\rho(w) = i - K$  proves that  $|Q| \leq 2$ . The case  $|Q| = 1$  can only occur if  $y_u^j, y_t^j$  and  $y_v^j$  are collinear for all  $j \leq p$ , an event that occurs with probability 0.

Proposition 5.4 is useful in order to show that certain configurations of nodes within  $\mathcal{T}$  can only occur with probability 0.

Example 5.5. Consider a subtree  $\mathcal{T}'$  of  $\mathcal{T}$  like the one in Fig. 6 embedded in  $\mathbb{R}^2$ , and suppose that all nodes at level  $u, w, t$  are marked  $\boxplus$ , and further that only one node within  $\alpha_1, \alpha_2$  is feasible, only one node within  $\alpha_3, \alpha_4$  is feasible, only one node within  $\alpha_7, \alpha_8$  is feasible, and  $\alpha_5, \alpha_6$  are both infeasible. This must be due to a distance  $d_{u'v}$  with  $u' \leq u$ . Consider now a circle  $C$  completely determined by its center at  $y_{u'}^1 = x(\alpha_1)_u$  and its radius  $d_{u'v}$ ; if  $C$  also contains the points at the nodes  $\alpha_1, \alpha_4, \alpha_8$  or the points at the nodes  $\alpha_2, \alpha_3, \alpha_7$  then we must have  $u' = u$ , in which case also one of  $\alpha_5, \alpha_6$  will be feasible (against the hypothesis). And the probability that  $C$  should contain the points at the nodes  $\alpha_1, \alpha_3, \alpha_8$  or  $\alpha_2, \alpha_4, \alpha_7$  is zero. Hence  $\mathcal{T}'$  can only occur with probability 0.

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

We now generalize Proposition 5.4 to vertices $u, v \in V$ with arbitrary rank difference.

Corollary 5.6. Consider a subtree $\mathcal{T}'$ of $\mathcal{T}$ consisting of $K + h + 1$ consecutive levels $i - K - h, \ldots, i$ (where $i \geq 2K + h$ and $h \geq 1$), rooted at a single node $\eta$ and such that all nodes at all levels are marked $\boxplus$. Let $p = 2^{K + h}$ and consider the set $Y' = \{y^j \mid j \leq p\}$ of partial realizations of $G$ at the leaf nodes of $\mathcal{T}'$. Let $u = \rho^{-1}(i - K - h)$ and $v = \rho^{-1}(i)$. Then with probability 1 there is a set $H^{uv} = \{r_j \mid j \leq 2^h\}$ of $2^h$ distinct positive reals such that $\| y'_u - y'_v \| \in H^{uv}$ for all $j \leq p$.

Proof. The proof of Proposition 5.4 can be generalized to span an arbitrary number of levels by induction on $q$. Two distances $r_{j_1}, r_{j_2} \in H^{uv}$ can only be equal by collinearity of some subsets of points, an event occurring with probability 0.

The next results show that, if $\{u, v\}$ is an edge of $E$ with rank difference greater than $K$, the distance $d_{uv}$ must belong to a certain finite set of values whenever the instance is a YES one.

Corollary 5.7. Let $y \in X$ and $v \in V \setminus V_0$ such that $\mathcal{T}(y, v)$ holds. If $\{u, w\} \in E$ with $u &lt; v &lt; w$ and $\rho(w) - \rho(u) &gt; K$ then $d_{uw} \in H^{uw}$ with probability 1.

Proof. Since $\mathcal{T}(y, v)$ holds, then the $^{\kappa}$ DMDGP instance is YES and there must exist at least two feasible nodes at level $\rho(w)$ in $\mathcal{T}$. If $d_{uw} \notin H^{uw}$ the probability for a given sphere to contain two arbitrary points in $\mathbb{R}^K$ is zero. Since the instance is a YES one, however, the BP algorithm does not prune all feasible nodes due to $d_{uw}$. By Corollary 5.6 the only remaining possibility (which therefore occurs with probability 1) is that $d_{uw} \in H^{uw}$.

Next, reflecting single points on realizations yields points that lie on the partial reflection of the whole realization.

Corollary 5.8. Let $y \in X$ and $v \in V \setminus V_0$ such that $\mathcal{T}(y, v)$ holds. If $u \in V$ with $u &gt; v$ then $R^v y_u$ belongs to a valid extension of $y[U_v]$.

Proof. If there is no edge $\{w, u\} \in E$ with $\rho(u) - \rho(w) &gt; K$ the result follows by Corollary 5.1. Otherwise, by Corollary 5.7, $d_{wu} \in H^{wu}$. As in the proof of Proposition 5.4, all pairs of points that are feasible with respect to $d_{wu}$ are reflections of each other with respect to $R^v$.

Finally, we state the main result of the section: if a $^{\kappa}$ DMDGP instance has a valid realization $y$ and $v$ is a vertex where a "valid branching" (in the sense of the $\mathcal{T}(y,v)$ assumption) takes place in the BP algorithm, then the partial reflection of $y$ with respect to $v$ is also a valid realization. This is surprising as the $\mathcal{T}(y,v)$ assumption only states that one of the branches at $v$ leads to $y$, whilst the other might end up at any other valid realization; Theorem 5.9 actually shows that the partial reflection of $y$ with respect to $v$ is valid.

Theorem 5.9. Let $y \in X$ and $v \in V \setminus V_0$ such that $\mathcal{T}(y, v)$ holds. Then $\bar{R}^v y \in X$ with probability 1.

Proof. We have to show that $\bar{R}^v y$ is a valid realization for $G$. Partition $E$ into three subsets $E_1, E_2, E_3$, where $E_1 = \{\{t, u\} \in E \mid t, u &lt; v\}$, $E_2 = \{\{t, u\} \in E \mid t, u \geq v\}$ and $E_3 = \{\{t, u\} \in E \mid t &lt; v \land u \geq v\}$. For $E_1$, by definition $\| (\bar{R}^v y)_t - (\bar{R}^v y)_u \| = \| ly_t - ly_u \| = \| y_t - y_u \| = d_{tu}$ as claimed. For $E_2$, $\| (\bar{R}^v y)_t - (\bar{R}^v y)_u \| = \| R^v y_t - R^v y_u \| = \| y_t - y_u \| = d_{tu}$ because $R^v$ is an isometry. For $E_3$, we aim to show that $\| ly_t - R^v y_u \| = d_{tu}$. Since $y \in X$, by Lemma 4.5 there is a feasible leaf node $\alpha$ with $x(\alpha) = y$. Because $\mathcal{T}(y, v)$, $\exists \eta \in \mathcal{V}_{\rho(v)-1}$ such that $x(\eta) = y[\gamma(v)]$ and $\mu(\beta) = \boxplus$ for all $\beta \in N^{+}(\eta)$. Let $\mathsf{p}(\alpha) \cap \mathcal{V}_{\rho(v)} = \{\beta\}$ for some $\beta \in N^{+}(\eta)$. Again by $\mathcal{T}(y, v)$, there is at least one feasible leaf node $\alpha'$ such that $\mathsf{p}(\alpha') \cap \mathcal{V}_{\rho(v)} = \{\beta'\}$ for some $\beta' \in N^{+}(\eta) \setminus \{\beta\}$. Let $\{\omega\} = \mathsf{p}(\alpha) \cap \mathcal{V}_{\rho(u)}$ and $\{\omega'\} = \mathsf{p}(\alpha') \cap \mathcal{V}_{\rho(u)}$. Because $\omega'$ is feasible, $\| x(\omega')_t - x(\omega')_u \| = d_{tu}$; because $\eta$ is an ancestor of both $\alpha$ and $\alpha'$ at level $\rho(v) - 1$ and $t &lt; v$, $\mathsf{p}(\alpha') \cap \mathcal{V}_{\rho(t)} = \mathsf{p}(\alpha) \cap \mathcal{V}_{\rho(t)}$, which implies that $x(\omega')_t = x(\omega)_t = y_t$. Thus, $\| y_t - y_u \| = d_{tu} = \| y_t - x(\omega')_u \|$. Furthermore, because $\beta' \in \mathsf{p}(\omega') \cap \mathcal{V}_{\rho(v)}, x(\omega')$ extends $x(\beta')$. By Algorithm 1, Steps 25 and 27, $\lambda(\beta) = 1 - \lambda(\beta')$. Because $\alpha$ is feasible, at every vertex $u' \in V$ such that $\rho(v) \leq \rho(u') &lt; \rho(u)$ the node $\theta \in \mathsf{p}(\alpha) \cap \mathcal{V}_{\rho(u')}$ has $f \in \{1, 2\}$ feasible subnodes; by Proposition 4.6, the node $\theta' \in \mathsf{p}(\alpha') \cap \mathcal{V}_{\rho(u')}$ also has $f$ feasible subnodes. If $f = 2$, by Corollary 5.8 and Lemma 5.3 it is possible to choose $\alpha'$ so that $\lambda(\theta') = 1 - \lambda(\theta)$ with probability 1; if $f = 1$ then by Algorithm 1, Steps 31 and 33, all feasible nodes inherit the same $\lambda$ value as their parents, so $\lambda(\theta') = 1 - \lambda(\theta)$. By Lemma 5.3, $x(\omega')_u = R^v y_u$ with probability 1. Hence $\| y_t - R^v y_u \| = d_{tu}$ as claimed.

## 6. Symmetry and the number of solutions

Our strategy for proving that YES instances of the $^{\kappa}$ DMDGP have power of two solutions with probability 1 is as follows. We map realizations $y \in X$ to binary sequences $\chi \in \{0, 1\}^n$ describing the "branching path" in the tree $\mathcal{T}$. We define a symmetry operation on $\chi$ by flipping its tail from a given component $i$ (this operation is akin to branching at level $i$, i.e. to a partial reflection applied to realizations). We show that the cardinality of the group of all such symmetries is a power of two by bijection with a set of binary sequences. Finally we prove that the cardinality of the symmetry group is the same as $|X|$.

For all leaf nodes $\alpha \in \mathcal{V}$ with $\mu(\alpha) = \boxplus$ let $\chi(\alpha) = (\lambda(\beta) \mid \beta \in \mathsf{p}(\alpha))$; since realizations in $X$ are also in correspondence with leaf $\boxplus$-nodes of $\mathcal{T}$ by Algorithm 1, Step 8, $\chi$ defines a relation on $X \times \{0, 1\}^n$.

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

Lemma 6.1. With probability 1, the relation $\chi$ is a function.

Proof. For $\chi$ to fail to be well-defined, there must exist a realization $x$ which is in relation with two distinct binary sequences $\chi', \chi''$, which corresponds to the discriminant of the quadratic equation in the proof of Lemma 4.1 taking value zero at some rank $&gt; K$, which happens with probability 0. $\square$

Let $\mathcal{Z} = \{\chi(y) \mid y \in X\}$. For $y \in X$ let $y^i$ be its subsequence $(x_1, \ldots, x_i)$. We extend $\chi$ to be defined on all such subsequences by simply setting $\chi^i = (\chi(y)_1, \ldots, \chi(y)_i)$; $\chi(y)$ is valid if $y$ is a valid realization.

Let $N = \{1, \dots, n\}$ and $g$ be the $n \times n$ binary matrix such that $g_{ij} = 1$ if $i \leq j$ and 0 otherwise (the upper triangular $n \times n$ all-1 matrix); let $g_i$ be its $i$-th row vector and $\Gamma = \{g_i \mid i \in N\}$. Consider the elementwise modulo-2 addition in the set $\mathbb{F}_2^n$ (denoted $\oplus$): this endows $\mathbb{F}_2^n$ with an additive group structure with identity $e = (0, \dots, 0)$ where each element is idempotent. Thus, $\mathcal{G} = (\mathbb{F}_2^n, \oplus) \cong C_2^n$. This group naturally acts on itself (and subsets thereof) using the same $\oplus$ operation. It is not difficult to prove that $\Gamma$ is a set of group generators for $\mathcal{G}$ and a linearly independent set of the vector space $\mathcal{V}$ given by $\mathcal{G}$ with scalar multiplication over $\mathbb{F}_2$. For all $S \subseteq N$, let

$$
g _ {S} = \bigoplus_ {i \in S} g _ {i},
$$

and define a mapping $\phi : \mathcal{P}(N) \to \mathcal{G}$ given by $\phi(S) = g_S$.

Lemma 6.2. $\phi$ is injective.

Proof. We show that for all $S, T \subseteq N$, if $g_S = g_T$ then $S = T$.

$$
\begin{array}{r c l}
\Rightarrow &amp; \bigoplus_ {i \in S} g _ {i} &amp; = \bigoplus_ {i \in T} g _ {i} \\
\Rightarrow &amp; \bigoplus_ {i \in S} g _ {i} \oplus \bigoplus_ {i \in T} g _ {i} ^ {- 1} &amp; = e
\end{array}
$$

$$
\begin{array}{l}
\text{idempotency} \quad \Rightarrow \quad \bigoplus_ {i \in S} g _ {i} \oplus \bigoplus_ {i \in T} g _ {i} = e \\
g _ {i} \oplus g _ {i} = g _ {i} ^ {2} \quad \Rightarrow \quad \bigoplus_ {i \in S \triangle T} g _ {i} \oplus \bigoplus_ {i \in S ^ {\prime}: T} g _ {i} ^ {2} = e \\
\text{idempotency} \quad \Rightarrow \quad \bigoplus_ {i \in S \triangle T} g _ {i} = e \\
\end{array}
$$

$$
\begin{array}{l}
\text{linear independence} \quad \Rightarrow \quad S \triangle T = \emptyset \\
\Rightarrow \quad S = T. \\
\end{array}
$$

This concludes the proof.

The following result shows essentially that groups of partial reflections have power of two cardinality.

Lemma 6.3. For all $H \subseteq \Gamma$, $|\langle H \rangle| = 2^{|H|}$.

Proof. The restriction of function $\phi$ to $\mathcal{P}(H)$ is injective by Lemma 6.2. Furthermore, each element $g$ of $\langle H\rangle$ can be written as $\bigoplus_{i\in S}g_i$ for some $S\subseteq H$ because $H$ is a spanning set for the vector space $H$ over $\mathbb{F}_2^n$, which is setwise equal to the group $\langle H\rangle$. Thus $\phi$ is surjective too. Hence $\phi$ is a bijection between $\mathcal{P}(H)$ and $\langle H\rangle$, which yields the result.

Let $I$ be the set of levels of $\mathcal{T}$ for which from all nodes with two feasible subnodes there is a path going to a feasible leaf through both. Let $L = \{g_i \in \Gamma \mid i \in I\}$ and $\Lambda = \langle L \rangle$ be the subgroup of $\mathcal{G}$ of "allowed partial reflections" generated by $L$. In the following (the main result of this section) we relate partial reflections to $\chi$ representations of valid realizations. We show that any valid realization, in its $\chi$ representation, generates the whole set of valid realizations by means of the action of the group of allowed partial reflections.

Theorem 6.4. If $\mathcal{Z} \neq \emptyset$, for all $\xi \in \mathcal{Z}$ we have $\xi \oplus \Lambda = \mathcal{Z}$ with probability 1.

Proof. $(\Rightarrow)$. We show that $\xi \oplus \Lambda \subseteq \mathcal{Z}$ with probability 1; because $\langle L\rangle = \Lambda$ it suffices to show that $\xi \oplus g_i\in \mathcal{Z}$ for an arbitrary $g_{i}\in L$, i.e. that there exists a valid realization $w\in X$ such that $\chi (w) = \xi \oplus g_{i}$. Let $y\in \chi^{-1}(\xi)$ and $v = \rho^{-1}(i)$ such that $\Upsilon (y,v)$, and define $w = \tilde{R}^{v}y$; by Theorem 5.9, $w\in X$. Let $\alpha^{\prime}$ be the leaf node of $\mathcal{T}$ such that $x(\alpha^{\prime}) = y$; by Lemma 4.5, there is a leaf node $\beta^\prime$ such that $x(\beta^{\prime}) = w$. We have to show that for all $\ell \geq i$ the node $\beta \in \mathsf{p}(\beta^{\prime})\cap \mathcal{V}_{\ell}$ is such that $\lambda (\beta) = 1 - \lambda (\alpha)$, where $\alpha$ is the node in $\mathsf{p}(\alpha^{\prime})\cap \mathcal{V}_{\ell}$. We proceed by induction on $\ell$. For $\ell = i$ this holds by Lemma 4.3. For $\ell &gt;i$, the induction hypothesis allows us to apply Lemma 5.3 and conclude that the event $\lambda (\alpha) = 1 - \lambda (\beta)$ occurs with probability 1.

$(\Leftarrow)$ Now we show that $\mathcal{Z} \subseteq \xi \oplus \Lambda$ with probability 1, i.e. for any $\eta \in \mathcal{Z}$ there is $g \in \Lambda$ with $\xi \oplus g = \eta$. We proceed by induction on $n$, which starts when $n = K + 1$: if $K + 1 \notin I$ then $|\mathcal{Z}| = 1$, $L = \emptyset$ and the theorem holds; if $K + 1 \in I$ then $|\mathcal{Z}| = 2$, $L = \{g_{K + 1}\}$ and the theorem holds. Now let $n &gt; K + 1$; for all $j \in \{K + 1, \dots, n - 1\}$ define $\mathcal{Z}^j = \{\xi^j \mid \xi \in \mathcal{Z}\}$

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

![img-7.jpeg](img-7.jpeg)
(a) Positions of the points on the plane.

![img-8.jpeg](img-8.jpeg)
(b) BP tree.
Fig. 7. The counterexample in the case  $K = 2$ . Realizations  $x_{5}^{(00)}, x_{5}^{(01)}$ , and  $x_{5}^{(11)}$  are valid, while  $x_{5}^{(10)}$  is not.

and  $L^i = \{g_\ell \in I^* \mid \ell \in I \land \ell \leq j\}$ . By the induction hypothesis, for all  $\xi' \in \mathcal{Z}^j$  ( $\xi' \oplus \langle L^i \rangle = \mathcal{Z}^j$ ). Now, either  $n \notin I$  or  $n \in I$ ; by Proposition 4.6, with probability 1 if  $n \notin I$  then nodes in  $\mathcal{V}_{n-1}$  can only have zero or one feasible subnode (let  $B_1^n$  be the set of all such feasible subnodes), and if  $n \in I$  then nodes in  $\mathcal{V}_{n-1}$  can only have zero or two feasible subnodes  $\beta$  (let  $B_2^n$  be the set of all such feasible subnodes). In the former case we let  $\mathcal{Z}^n = \{\xi(x(\beta)) \mid \beta \in B_1^n\}$  and  $L^n = L^{n-1}$ ; in the latter we let  $\mathcal{Z}^n = \{\xi(x(\beta)) \mid \beta \in B_2^n\}$  and  $L^n = L^{n-1} \cup \{g_n\}$ . In both cases it is easy to verify that the theorem holds for  $\mathcal{Z}^n$ ,  $L^n$ : in the former case it follows by the induction hypothesis, and in the latter case it follows because  $g_n = (0, \dots, 0, 1)$ , namely, if  $\eta \in \mathcal{Z}$  and  $n \in I$  then take  $\xi = \eta \oplus g_n$  (the result follows by idempotency of  $g_n$ ).

The main result of the paper is now simply a corollary of Theorem 6.4.

Corollary 6.5. If a  $^{\kappa}$  DMDGP instance is YES,  $|X|$  is a power of two with probability 1.

Proof. By Lemma 6.1  $\chi$  is a function with probability 1. Let  $x, x' \in X$  be distinct; then by Algorithm 1, Steps 25, 27, 31, and 33, the map  $\chi : X \to \mathcal{Z}$  is injective. By definition of  $\mathcal{Z}$  it is also surjective, hence  $|X| = |\mathcal{Z}|$ . By Theorem 6.4  $|\mathcal{Z}| = |\chi \oplus \Lambda|$  for all  $\chi \in \mathcal{Z}$  with probability 1. It is easy to show that  $|\chi \oplus \Lambda| = |\Lambda|$ , so by Lemma 6.3  $|X|$  is a power of two with probability 1.

# 7. Counterexamples

# 7.1. Disproving the "power of two" conjecture

We first discuss a class of counterexamples to the conjecture that all  $^{\kappa}$  DMDGP instances have a number of solutions which is a power of two (also see Lemma 5.1 in [21]). All these counterexamples are hand-crafted and have the property that two distinct realizations  $x, x'$  have at least a level  $i$  where  $x_i = x'_i$ , which is an event which happens with probability 0. For any  $K \geq 1$ , let  $n = K + 3$ ,  $V = \{1, \ldots, n\}$ ,  $E = \{\{i, j\} \mid 0 &lt; i - j \leq K\} \cup \{\{1, n\}\}$  and  $d_{\bar{q}} = 1$  for all  $\{i, j\} \in E$ . The first  $n - 2 = K + 1$  vertices can be realized as the vertices of a regular simplex in dimension  $K$ ; then either  $x_{n-1} = x_1$  or  $x_{n-1}$  is the symmetric position from  $x_1$  with respect to the hyperplane through  $\{x_2, \ldots, x_{n-2}\}$ . In the first case, the two positions for  $x_n$  are valid, in the second only  $x_n = x_2$  is possible (see Fig. 7 for the 2-dimensional case), yielding a YES instance where  $|X| = 6$ .

# 7.2. Necessity of immediate predecessors

Lastly, Fig. 8 shows an example where the  $\langle \mathrm{ii} \rangle \Rightarrow \langle \mathrm{i} \rangle$  implication of Lemma 5.3 fails for instances in DDGP  $\setminus$  DMDGP. This shows that any generalization of our result to the DDGP must be nontrivial. Let  $V = \{1, \ldots, 6\}$  (the graph drawing is the same as the realization in  $\mathbb{R}^2$ ). The points  $5', 6'$  linked with dashed lines show alternative placements for the corresponding vertices. Let  $U_5 = \{3, 4\}$  and  $U_6 = \{1, 2\}$ . The line through the points 3, 4 does not provide a valid reflection mapping 6 to  $6'$ . This happens because  $U_6$  does not consist of the two immediate predecessors of vertex 6.

# 8. Application of the theory to protein conformation

The "interval BP algorithm" (iBP) is the adaptation of the BP algorithm to a more realistic setting for protein backbones [23]. For  $v \in V$ , we shall denote the vertex of rank  $\rho(v) - h$  (for some  $h &lt; \rho(v)$ ) by  $v - h$  (in other words, we identify the vertex labels with their own ranks in the order). We note that throughout this section  $K$  is fixed to 3; we still use the symbolic notation  $K$  to emphasize results that hold for whatever value of  $K$ .

As mentioned in the introduction, covalent bonds and angles are known reasonably precisely. These allow us to consider the distances  $d_{v - 2,v}$  and  $d_{v - 1,v}$  as precise for any  $v$  with rank greater than  $K$ . The data provided by NMR actually consists of a

---

L. Liberti et al./ Discrete Applied Mathematics 165 (2014) 213-232

![img-9.jpeg](img-9.jpeg)
Fig. 8. A counterexample to Lemma 5.3 applied to DDGP  $\searrow$  DMDGP.

frequency measurement for a triplet  $(a, b, \delta)$  where  $a, b$  are atom types (e.g. H, C, O and so on) and  $\delta$  is a distance value [15]. In other words, NMR gives experimental evidence that the distance  $\delta$  is to be expected between atoms of type  $a, b$  with a certain probability. NMR specialists then solve an assignment problem so that they can decide which pair(s) of atoms of type  $a, b$  should be assigned the distance  $\delta$ . This procedure gives rise to errors, due to which all distances are actually reported as intervals rather than precise values. We emphasize two important observations:

1. By default, NMR experiments are rigged up to find the distances between pairs of hydrogen atoms only [43]. Distances between pairs of atoms of different types are possible but they require more work, and are prone to more errors.
2. The NMR machinery has a resolution limit of around about one decimal digit (in A units); if the distance interval is on a similar scale, only a finite number of values belonging to the intervals can actually be measured [39].

We initially exploited Observation 1 by defining certain vertex orders that only include hydrogens [27]; following a similar methodology, we later defined another vertex order (called the iBP order) which also takes into account Observation 2 [23]. More precisely, for each  $v \in V$  of rank greater than 3, this order is such that:

-  $d_{v - 2,v}, d_{v - 1,v}$  are single values in  $\mathbb{R}_+$ ;
-  $d_{v - 3,v}$  is either a single value, or a finite set  $\{d_{uv}^1,\dots ,d_{uv}^q\}$  (with  $q &gt; 1$ ) of values in  $\mathbb{R}_+$ ,
-  $d_{uv}$  is either a single value, or an interval  $[d_{uv}^k, d_{uv}^O]$  for all  $\rho(u) &lt; \rho(v) - 3$ .

Essentially, this order is such that all interval distances are only used for pruning purposes, whilst the discretization of the search space occurs differently because  $d_{v - 3,v}$  could be a finite set of values.

# 8.1. Discretization with a finite set of values

In this section we go through the previous results and adapt them to the iBP setting insofar as  $d_{v-3,v}$  might be a set instead of a single value. We still consider precise distances for pruning rather than intervals (this further extension will be tackled in Section 8.2). Before going into the details, the intuitive explanation as to why most results still hold (albeit with some modifications) is that one can see the iBP search tree as the union of  $q^k$  different BP trees, where  $k$  is the number of vertices  $v \in V$  such that  $d_{v-3,v}$  is a set. Thus the iBP tree inherits many of the properties of all the BP trees that compose it.

# 8.1.1. Changes to Section 4

Assume  $v \in V$  is such that  $d_{v-3,v}$  is a set of  $q$  values and every other distance is precise. Then Lemma 4.1 changes as follows.

Lemma 8.1. If  $|N(v) \cap \gamma(v)| = 3$  then there are at most  $2q$  distinct extensions of  $x$  that are valid for  $G^v$ . If one valid extension exists, then with probability 1 there are exactly  $2q$  distinct valid extensions.

Proof. We shall only point out the differences with the proof of Lemma 4.1. We choose the vertex  $w \in U_v$  such that  $w - 3$ ; then Eq. (5) becomes  $\forall u \in U_v \setminus \{w\} (\| z - x_u \| = d_{uv}) \land (\| z - x_w \| \in d_{uv})$ , so that the second line in Eq. (6) becomes  $\| z \|^2 - 2x_w \cdot z + \| x_w \|^2 \in d_{uv}^2$ , where  $d_{uv}^2$  denotes the set  $\{(d_{uv}^j)^2 \mid j \leq q\}$ . We rewrite this as  $q$  different quadratic equations

$$
\| z \| ^ {2} - 2 x _ {w} \cdot z + \| x _ {w} \| ^ {2} = d _ {w v} ^ {j} \tag {10}
$$

for each  $j \leq q$ . From each of these, by Lemma 4.1 we obtain two distinct values for  $z$  with probability 1, which concludes the proof.

Lemma 4.2 involves adjacent but not immediate predecessors. The statement is unchanged, but the proof is slightly different.

Lemma 8.2. If  $|N(v) \cap \gamma(v)| &gt; 3$  then, with probability 1, there is at most 1 extension of  $x$ .

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

Proof. Since $U_v \subseteq N(v) \cap \gamma(v)$, only considering vertices in $U_v$ we are in the situation of Lemma 8.1 and we have $2q$ possible extensions of $x$. Now consider an adjacent predecessor $u$ of $v$ which is not in $U_v$: with probability 1, the distances $\| x_u - z^{1j} \|, \| x_u - z^{1j} \|$ are all distinct (for all $j \leq q$). Hence, at most one can be equal to $d_{uv}$.

Lemma 4.3 can be adapted in the following way.

Lemma 8.3. Let $\tilde{x}$ be any valid realization of $G[U_v]$ and $\{z^{1j}, z^{2j} \mid j \leq q\} \subseteq \mathbb{R}^3$ be the $2q$ positions for vertex $v$ given by Lemma 8.1, with $z^{1j}, z^{2j}$ arising from the Eq. (10) indexed by $j$. Then $z^{2j}$ is a reflection of $z^{1j}$ with respect to the hyperplane through the 3 points of $\tilde{x}$.

The proof is the same as that of Lemma 4.3, for each $j \leq q$.

The change to these lemmata impacts Algorithm 1 in several respects. The most important change is that whenever $d_{v-3,v}$ is a set and $v$ has exactly 3 adjacent predecessors, the BP tree node $\alpha$ at level $v-1$ has exactly $2q$ subnodes $\alpha^{ij}$ (for $i \in \{1,2\}$ and $j \leq q$) at level $v$. Accordingly, for all $i \in \{1,2\}$ and $j \leq q$, $\lambda(\alpha^{ij})$ is defined as follows:

$$
\lambda \left(\alpha^{ij}\right) = \begin{cases}
-j &amp; \text{if } a_v^\top z^{ij} \leq a_{v0} \\
j &amp; \text{if } a_v^\top z^{ij} &gt; a_{v0}.
\end{cases} \tag{11}
$$

The set $Z$ in Step 17 of Algorithm 1 is $\{z^{ij} \mid i \in \{1,2\} \land j \leq q\}$. The test at Step 23 is $|Z| = 2q$. The body of the corresponding if changes in the following way:

for $j \leq q$ do

if $a_v^\top z^{1j} \leq a_{v0}$ then

Let $\alpha^{1j} = (z^{1j}, -j, \boxplus), \alpha^{2j} = (z^{2j}, j, \boxplus)$

else

Let $\alpha^{1j} = (z^{1j}, j, \boxplus), \alpha^{2j} = (z^{2j}, -j, \boxplus)$

end if

end for

It is easy to show that the full $i$ BP search tree is no longer a binary one: every time $d_{v-3,v}$ is a set, level each node at level $v-1$ has $2q$ subnodes. The total number of $i$ BP nodes depends on the incidence of set distances within the $^{\mathrm{K}}$ DMDGP vertex order.

Lemma 4.5 does not change. Proposition 4.6 changes as follows.

Proposition 8.4. With probability 1, there is no level $i \leq n$ having two distinct feasible nodes $\beta, \theta \in \mathcal{V}_i$ such that $|\{\alpha \in N^+(\beta) \mid \mu(\alpha) = \boxplus\}| = 1$ and $|\{\alpha \in N^+(\theta) \mid \mu(\alpha) = \boxplus\}| = 2q$.

The proof is trivially adapted from that of Proposition 4.6.

## 8.1.2. Changes to Section 5

The statement of Corollary 5.1 changes as follows.

Corollary 8.5. Let $\alpha \in \mathcal{V}_{i-1}$ for some $i &gt; 1$, $v = \rho^{-1}(i)$ and $N^{+}(\alpha) = \{\eta^{ij} \mid i \in \{1,2\} \land j \leq q\}$ with $\mu(\eta^{ij}) = \boxplus$ for all $i \in \{1,2\}$, $j \leq q$. Then $x(\eta^{1j})_v = R^v x(\eta^{2j})_v$.

Only Condition (ii) changes in Lemma 5.3:

(ii) $\forall i \leq j \leq \ell$ if $d_{uw}$ is a single value, then Condition (ii) in Lemma 5.3 holds; if $d_{uw}$ is a set, $\lambda(\eta^v) = -\lambda(\beta^v)$.

The proof is an easy adaptation of that of Lemma 5.3: every time $\lambda$ is mentioned, it suffices to verify whether the corresponding distance $d_{uw}$ is a value or a set, and use the correct definition (either Eq. (9) or Eq. (11)). Proposition 5.4 changes its statement as follows.

Proposition 8.6. Consider a subtree $\mathcal{T}'$ of $\mathcal{T}$ consisting of $K + 2$ consecutive levels $i - K - 1, \ldots, i$ (where $i \geq 2K + 1$), rooted at a single node $\eta$ and such that all nodes at all levels are marked $\boxplus$. Let $u = \rho^{-1}(i - K - 1)$, $w = \rho^{-1}(i - K)$, $t = \rho^{-1}(i - 1)$ and $v = \rho^{-1}(i)$. If $p$ is the number of leaf nodes of $\mathcal{T}'$, let $Y' = \{y^j \mid j \leq p\}$ be the set of partial realizations of $G$ at the such leaf nodes.

1. If $d_{ut}$ and $d_{wv}$ are single values, Proposition 5.4 holds.
2. If $d_{ut}$ is a set of $q$ values and $d_{wv}$ is a single value or vice versa, let $p = 2^{K + 1}q$; with probability 1 there is a set $H^{uv} \subseteq \mathbb{R}_+$ with $|H^{uv}| = 2q$ such that $\| y_u^j - y_v^j \| \in H^{uv}$ for each $j \leq p$.
3. If both $d_{ut}$, $d_{wv}$ are sets of $q$ values, let $p = 2^{K + 1}q^2$; with probability 1 there is a set $H^{uv} \subseteq \mathbb{R}_+$ with $|H^{uv}| = 2q^2$ such that $\| y_u^j - y_v^j \| \in H^{uv}$ for each $j \leq p$.

---

L. Liberti et al./ Discrete Applied Mathematics 165 (2014) 213-232

![img-10.jpeg](img-10.jpeg)
Fig. 9. Proof of Proposition 5.4 adapted to the iBP setting in  $\mathbb{R}^2$ , with  $d_{ut}$  having cardinality 2 and  $d_{uv}$  being a single value. The arcs of circles  $r_k^j (j \leq 2, h \in \{1, 2\})$  centered at  $u$  mark the four possible distances between the position of  $u$  and that of  $v$ .

Proof. The situation in  $\mathbb{R}^2$  is shown in Fig. 9 for the case where  $d_{ut}$  is a set of 2 values and  $d_{uv}$  is a single value. Only two edges, namely  $\{u,t\}$  and  $\{w,v\}$ , can be weighted by sets of distance values in the subtree  $\mathcal{T}'$ . If  $\{u,t\}$  is weighted by a set but  $\{w,v\}$  is not, then at level  $t$  we have  $2q$  as many nodes as in the previous level; if  $\{w,v\}$  is weighted by a set but  $\{u,t\}$  is not, then at level  $v$  we have  $2q$  as many nodes as in level  $t$ ; if both are weighted by sets, then at both levels the increase in the number of nodes is  $2q$ -fold. At any level  $s\in \{t,v\}$  where there is an increase of  $2q$  nodes w.r.t. the previous level, each pair of nodes  $\alpha^{1j},\alpha^{2j}$  with  $\lambda (\alpha^{1j}) = -\lambda (\alpha^{2j})$  (for  $j\leq p$ ) is such that  $x(\alpha^{1j})_s$  and  $x(\alpha^{2j})_s$  are reflections through  $R^w$  by the adaptation of Lemma 5.3 (above). Hence each such pair yields a subtree of  $\mathcal{T}'$  which, by Proposition 5.4, with probability 1 allows  $\| y_u^j -y_v^j\|$  to only take values in a set  $S\subseteq \mathbb{R}_{+}$  such that  $|S| = 2$ . Observe that any pair of such sets  $S$  will have non-empty intersection with probability 0. The result follows.

Corollary 5.6 can be extended to the iBP setting by remarking that the set  $H^{uv}$  (appearing in its statement) has cardinality  $2^h q^k$ , where  $k$  is the number of vertices  $t$  between  $u$  and  $v$  such that  $d_{t-3,t}$  is a set. Corollaries 5.7 and 5.8 are unchanged.

The statement of Theorem 5.9 is unchanged, but the proof needs to be adapted. The statement  $\lambda (\beta) = 1 - \lambda (\beta^{\prime})$  must be changed to  $\lambda (\beta) = -\lambda (\beta^{\prime})$  by the adaptation of Algorithm 1 to the iBP given above in this section. The statement  $f\in \{1,2\}$  must be changed to "  $f\in \{1,2\}$  if  $d_{u^{\prime} - 3,u^{\prime}}$  is a single value, and  $f\in \{1,2q\}$  if  $d_{u^{\prime} - 3,u^{\prime}}$  is a set". The motivation for  $\theta^\prime$  to also have  $f$  subnodes is by Proposition 4.6 if  $d_{u^{\prime} - 3,u^{\prime}}$  is a single value, and by Proposition 8.4 if  $d_{u^{\prime} - 3,u^{\prime}}$  is a set. Between the cases  $f = 2$  and  $f = 1$ , the case  $f = 2q$  must be added, with the comment that it is possible to choose  $\alpha^\prime$  so that  $\lambda (\theta^{\prime}) = -\lambda (\theta)$  with probability 1 (by Lemma 5.3 adapted to the iBP setting). The case  $f = 1$  must be changed so that if  $d_{u^{\prime} - 3,u^{\prime}}$  is a set, then  $\lambda (\theta^{\prime}) = -\lambda (\theta)$ .

# 8.1.3. The number of solutions

Since the iBP branches over  $2q$  possibilities at some levels, in general the number of possible realization is no longer a power of two, but rather a number  $2^{\ell}q^{k}$  for some integers  $\ell, k$ . Because of this fundamental difference, we shall not attempt to adapt each result in Section 6 to the iBP setting, but rather propose a somewhat different development.

For all  $v \in V$ , let  $\psi_v = \{0\}$  whenever  $v \leq K$ ,  $\psi_v = \{0, 1\}$  whenever  $d_{v-3,v}$  is a single value, and  $\psi_v = \{0, \dots, 2q-1\}$  whenever  $d_{v-3,v}$  is a set. Accordingly, we define  $\chi$  to map realizations in  $X$  to sequences indexed by  $V$ , in such a way that  $\chi_v \in \psi_v$  for all  $v \in V$ . More precisely, if  $\alpha$  is an iBP tree node such that  $x(\alpha)_v = y_v$  and  $d_{v-3,v}$  is a set, then  $\chi_v = \lambda(\alpha) - 1$  whenever  $\lambda(\alpha) &gt; 0$ , and  $\chi_v = q - \lambda(\alpha) - 1$  whenever  $\lambda(\alpha) &lt; 0$ . With this definition, the set  $\Xi = \{\chi(y) | y \in X\}$  can, with probability 1, be endowed with an Abelian group structure:  $\chi^\dagger \oplus \chi^2 = (\chi_v^\dagger \oplus \chi_v^2 | v \in V)$ , where  $\oplus$  denotes addition modulo  $|\psi_v|$ . We remark that the group  $\mathcal{G} = \langle \Xi, \oplus \rangle$  is a subgroup of  $\tilde{\mathcal{G}} = \prod_{v \in V} C_{|\psi_v|}$ . Furthermore, since the underlying set of  $\mathcal{G}$  is  $\Xi, \mathcal{G}$  acts on itself: this yields a regular action, which is transitive by definition. Thus, for any  $\chi \in \Xi$ , we have  $\chi \oplus \Xi = \Xi$ .

For all  $v \in V$  such that  $\rho(v) &gt; 3$ , if  $d_{v-3,v}$  is a single value, let  $\iota_v = 1$ ; if  $d_{v-3,v}$  is a set, let  $\iota_v = q$ . Consider the set  $\tilde{\Lambda}$  of the following elements of  $\tilde{\mathcal{G}}$ :

-  $\pi_v = (0, \dots, 0, \iota_v, \dots, \iota_n)$  for all  $v$  such that  $\rho(v) &gt; 3$ ;
-  $\sigma_v = (0, \dots, 0, 1_v, 0, \dots, 0)$ , where the 1 is in the  $\rho(v)$ -th position, for all  $v$  such that  $d_{v-3,v}$  is a set.

It is easy to show that  $\tilde{\Lambda}$  generates the whole of  $\tilde{\mathcal{G}}$ . Now, similarly to Section 6, we find a subset  $\Lambda$  of  $\tilde{\Lambda}$  which generates  $\mathcal{G}$ . More precisely:

- if  $v \in V$  is such that  $\rho(v) &gt; 3$ ,  $d_{v - 3,v}$  is a single value, and for all iBP nodes at level  $v$  with two feasible subnodes, both subnodes are on paths continuing to leaf nodes, then  $\pi_v \in \Lambda$ ;

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

![img-11.jpeg](img-11.jpeg)
Fig. 10. The interval distance $d_{uv}$ makes 3 out of 4 possible values from $H^{uv}$ feasible, namely $r_1^1, r_1^2, r_2^2$ (and hence 12 out of 16 possible positions for vertex $v$).

- if $v \in V$ is such that $\rho(v) &gt; 3$, $d_{v - 3, v}$ is a set and for all iBP nodes at level $v$ with $2q$ feasible subnodes, both subnodes are on paths continuing to leaf nodes, then $\pi_v, \sigma_v \in \Lambda$.

We are now in a position to prove a theorem whose significance is similar to that of Theorem 6.4.

Theorem 8.7. With probability 1, $\langle \Lambda ,\oplus \rangle = \tilde{y}$

Proof. That $\langle \Lambda, \oplus \rangle \subseteq \tilde{y}$ follows by definition: $\pi_v$ encodes a valid branching, and the elements $\pi_v \oplus (\sigma_v \oplus \dots \sigma_v)$ represent the $2q$ branches in the case of set distances. The other direction is similar to the $\langle \Leftarrow \rangle$ direction of the proof of Theorem 6.4. Replace "zero or two feasible subnodes" by "zero, two or $2q$ feasible subnodes", and $\{g_n\}$ by $\{\pi_n\}$ if there are two feasible subnodes at level $n$, or by $\{(\pi_n)^k \oplus (\sigma_v)^h \mid k \in \{0,1\} \land 0 \leq h \leq q - 1\}$ otherwise, where $\alpha^f = \sum_{s &lt; f} \alpha$ for any integer $f$ and $\alpha \in \tilde{y}$.

Theorem 8.7, like Theorem 6.4, can be used to compute the whole of $X$ having just one $x \in X$ and the partial reflection operators in $\tilde{y}$.

Let $M$ be the number of vertices $v$ with $\rho(v) &gt; 3$ such that $d_{v-3,v}$ is a set. It is easy to see that the cardinality of $\langle \Lambda, \oplus \rangle$ is $2^{\ell}q^{k}$ for some $\ell \leq n-3$ and $k \leq M$, since it is a subgroup of $\tilde{\mathcal{G}}$. By Theorem 8.7, with probability 1, that is the same cardinality of $\tilde{\mathcal{G}}$, and therefore of $\Xi$.

## 8.2. Pruning by intervals

If $z \in \mathbb{R}^K$ is a candidate position for vertex $v$, and $u$ is a non-immediate adjacent predecessor of $v$ (i.e., such that $\rho(u) &lt; \rho(v) - 3$), then $d_{uv}$ may be the interval $[d_{uv}^L, d_{uv}^G]$, and the test on Line 19 of Algorithm 1 becomes $\| x(\beta)_u - z \| \notin [d_{uv}^L, d_{uv}^G]$.

It was shown in Proposition 5.4, Corollary 5.6, Proposition 8.6, and the adaptation of Corollary 5.6 to distance sets in Section 8.1.2 that with probability 1, for $u, v \in V$ with $\rho(v) - \rho(u) &gt; K$, there are sets $H^{uv} \subseteq \mathbb{R}_+$ of an even number of possible distance values that $\| y_u - y_v \|$ can take for any $y \in X$. Furthermore, with probability 1, for any distance value $r \in H^{uv}$ and partial realization $(y_1, \ldots, y_u)$ there are two possible extensions $y, y'$ to a partial realization on $G[\gamma(v) \cup \{v\}]$ such that $\| y_u - y_v \| = r$; moreover, it turns out that we have $y_w = R^{u+1} y_w'$ for all $w \in V$ such that $\rho(u) &lt; \rho(w) \leq \rho(v)$. By the results of Section 5, a single value $d_{uv} = r$ weighting an edge $\{u, v\} \in E$ prunes all possible realizations for $v$ apart from those two for which $\| y_u - y_v \| = r$. Now, if $d_{uv}$ is an interval, then it can contain any number of consecutive values in $H^{uv}$, and as such it can prune any even number of possible realizations for $v$ (including zero), as shown graphically in Fig. 10.

It is therefore difficult to envisage an aprioristic calculation of the number of solutions $|X|$ of a given instance in this case. For YES instances, it is clear that the result obtained in Section 8.1.3 provides a lower bound to this number; and evidently the cardinality of $\tilde{\mathcal{G}}$ provides an upper bound. Thus, there exist $\ell \leq n - 3$ and $k \leq M$ (which depend on vertices $v$ adjacent to edges $\{u, v\}$ where $u &lt; v - 3$) such that:

$$
2^{\ell} q^{k} \leq |X| \leq 2^{n-3} q^{M}. \tag{12}
$$

## 8.3. Distance sets of different cardinalities

Although we assumed that, for all $v \in V$ such that $\rho(v) &gt; 3$ and $d_{v-3,v}$ is a set of values, the cardinality of $|d_{v-3,v}|$ is always $q$, this may fail to be the case; we shall let $q_v = |d_{v-3,v}|$ for all such $v$ (for the other $v \in V$, we define $q_v = 1$). Furthermore, we let $\epsilon \in \{0,1\}^n$ such that $\epsilon_v = 0$ for all $\rho(v) \leq 3$ and all $\rho(v)$ such that $d_{v-3,v}$ is a single value. The other

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

components of $v$ will depend on whether $v$ has $K$ or more than $K$ adjacent predecessors. It is not difficult to adapt the results of Section 8.1 in order to show that Eq. (12) becomes:

$$
2^{\ell} \prod_{v \in V} q_{v}^{\epsilon_{v}} \leq |X| \leq 2^{n-3} \prod_{v \in V} q_{v}. \tag{13}
$$

Both lower and upper bounds are modifications of Eq. (12), where $q^k$ has been replaced by the product of the values $q_v$ when $v$ ranges over the relevant vertices of $V$.

## 9. Conclusion

This paper presents a mathematical theory of non-unique rigidity in connection to symmetry of certain graphs arising in protein conformation problems from NMR data. Further work will focus on a computational application of these results, where we shall exploit Theorem 6.4 in order to compute the set of all incongruent realizations of a given $^{\mathrm{K}}$ DMDGP graph using just one realization thereof.

## References

[1] A. Barvinok, Problems of distance geometry and convex properties of quadratic maps, Discrete and Computational Geometry 13 (1995) 189-202.
[2] B. Berger, J. Kleinberg, T. Leighton, Reconstructing a three-dimensional model with arbitrary errors, Journal of the ACM 46 (2) (1999) 212-235.
[3] P. Biswas, T. Lian, T. Wang, Y. Ye, Semidefinite programming based algorithms for sensor network localization, ACM Transactions on Sensor Networks 2 (2006) 188-220.
[4] L. Blumenthal, Theory and Applications of Distance Geometry, Oxford University Press, Oxford, 1953.
[5] T. Brady, C. Watt, On products of Euclidean reflections, American Mathematical Monthly 113 (2006) 826-829.
[6] R. Connelly, Generic global rigidity, Discrete &amp; Computational Geometry 33 (2005) 549-563.
[7] I. Coope, Reliable computation of the points of intersection of $n$ spheres in $\mathbb{R}^n$, The Australian and New Zealand Industrial and Applied Mathematics Journal 42 (2000) C461-C477.
[8] L. Cremona, Le Figure Reciproche Nella Statica Grafica, G. Bernardoni, Milano, 1872.
[9] G. Crippen, T. Havel, Distance Geometry and Molecular Conformation, Wiley, New York, 1988.
[10] J. Dattorro, Convex Optimization and Euclidean Distance Geometry, Meboo, Palo Alto, 2005.
[11] Q. Dong, Z. Wu, A geometric build-up algorithm for solving the molecular distance geometry problem with sparse distance data, Journal of Global Optimization 26 (2003) 321-333.
[12] T. Eren, D. Goldenberg, W. Whiteley, Y. Yang, A. Morse, B. Anderson, P. Belhumeur, Rigidity, computation, and randomization in network localization, in: IEEE Infocom Proceedings, 2004, pp. 2673-2684.
[13] J. Graver, Rigidity matroids, SIAM Journal on Discrete Mathematics 4 (1991) 355-368.
[14] J. Graver, B. Servatius, H. Servatius, Combinatorial Rigidity, American Mathematical Society, 1993.
[15] H. Gunther, NMR Spectroscopy: Basic Principles, Concepts, and Applications in Chemistry, Wiley, New York, 1995.
[16] B. Hendrickson, Conditions for unique graph realizations, SIAM Journal on Computing 21 (1) (1992) 65-84.
[17] L. Henneberg, Die Graphische Statik der Starren Systeme, Teubner, Leipzig, 1911.
[18] A.L.-S. John, Geometric constraint systems with applications in cad and biology, Ph.D. Thesis, University of Massachusetts at Amherst, 2008.
[19] R. Kang, T. Müller, Sphere and dot product representations of graphs, in: Proceedings of SCG11, ACM, 2011, pp. 308-314.
[20] C. Lavor, J. Lee, A.L.-S. John, L. Liberti, A. Mucherino, M. Sviridenko, Discretization orders for distance geometry problems, Optimization Letters 6 (2012) 783-796.
[21] C. Lavor, L. Liberti, N. Maculan, The discretizable molecular distance geometry problem, Tech. Rep. q-bio/0608012, arXiv, 2006.
[22] C. Lavor, L. Liberti, N. Maculan, A. Mucherino, The discretizable molecular distance geometry problem, Computational Optimization and Applications 52 (2012) 115-146.
[23] C. Lavor, L. Liberti, A. Mucherino, The interval Branch-and-Prune algorithm for the discretizable molecular distance geometry problem with inexact distances, Journal of Global Optimization (2011). http://dx.doi.org/10.1007/s10898-011-9799-6.
[24] C. Lavor, A. Mucherino, L. Liberti, N. Maculan, An artificial backbone of hydrogens for finding the conformation of protein molecules, in: Proceedings of the Computational Structural Bioinformatics Workshop, IEEE, Washington, DC, USA, 2009, pp. 152-155.
[25] C. Lavor, A. Mucherino, L. Liberti, N. Maculan, Computing artificial backbones of hydrogen atoms in order to discover protein backbones, in: Proceedings of the International Multiconference on Computer Science and Information Technology, IEEE, Mragowo, Poland, 2009, pp. 751-756.
[26] C. Lavor, A. Mucherino, L. Liberti, N. Maculan, Discrete approaches for solving molecular distance geometry problems using NMR data, International Journal of Computational Bioscience 1 (2010) 88-94.
[27] C. Lavor, A. Mucherino, L. Liberti, N. Maculan, On the computation of protein backbones by using artificial backbones of hydrogens, Journal of Global Optimization 50 (2011) 329-344.
[28] L. Liberti, C. Lavor, N. Maculan, A branch-and-prune algorithm for the molecular distance geometry problem, International Transactions in Operational Research 15 (2008) 1-17.
[29] L. Liberti, C. Lavor, A. Mucherino, N. Maculan, Molecular distance geometry methods: from continuous to discrete, International Transactions in Operational Research 18 (2010) 33-51.
[30] L. Liberti, B. Masson, J. Lee, C. Lavor, A. Mucherino, On the number of solutions of the discretizable molecular distance geometry problem, in: Combinatorial Optimization, Constraints and Applications, COCOA11, in: LNCS, vol. 6831, Springer, New York, 2011, pp. 322-342.
[31] K. Menger, Untersuchungen über allgemeine Metrik, Mathematische Annalen 100 (1928) 75-163.
[32] A. Mucherino, C. Lavor, The branch and prune algorithm for the molecular distance geometry problem with inexact distances, in: Proceedings of the International Conference on Computational Biology, in: World Academy of Science, Engineering and Technology, vol. 58, 2009, pp. 349-353.
[33] A. Mucherino, C. Lavor, L. Liberti, Exploiting symmetry properties of the discretizable molecular distance geometry problem, Journal of Bioinformatics and Computational Biology 10 (2012) 1242009 (1-15).
[34] A. Mucherino, C. Lavor, L. Liberti, The discretizable distance geometry problem, Optimization Letters 6 (2012) 1671-1686.
[35] A. Mucherino, C. Lavor, L. Liberti, N. Maculan, On the definition of artificial backbones for the discretizable molecular distance geometry problem, Mathematica Balkanica 23 (2009) 289-302.
[36] A. Mucherino, C. Lavor, L. Liberti, E.-G. Talbi, A parallel version of the branch &amp; prune algorithm for the molecular distance geometry problem, in: ACS/IEEE International Conference on Computer Systems and Applications, AICCSA10, IEEE, Hammamet, Tunisia, 2010, pp. 1-6.
[37] A. Mucherino, C. Lavor, L. Liberti, E.-G. Talbi, On suitable parallel implementations of the branch &amp; prune algorithm for distance geometry, in: Proceedings of the Grid5000, Spring School, Lille, France, 2010.

---

L. Liberti et al. / Discrete Applied Mathematics 165 (2014) 213-232

[38] A. Mucherino, L. Liberti, C. Lavor, MD-jeep: an implementation of a branch-and-prune algorithm for distance geometry problems, in: K. Fukuda, J. van der Hoeven, M. Joswig, N. Takayama (Eds.), Mathematical Software, in: LNCS, vol. 6327, Springer, New York, 2010, pp. 186–197.

[39] M. Nilges, M. Macias, S. O'Donoghue, H. Oschkinat, Automated NOESY interpretation with ambiguous distance restraints: the refined NMR solution structure of the Pleckstrin homology domain from $\beta$-spectrin, Journal of Molecular Biology 269 (1997) 408–422.

[40] C. Saviotti, Nouvelles méthodes pour le calcul des travaux réticulaires, in: Appendix to L. Cremona, Les Figures Réciproques en Statistique Graphique, Gauthier-Villars, Paris, 1885, pp. 37–100.

[41] C. Saviotti, La Statica Grafica: Lezioni, U. Hoepli, Milano, 1888.

[42] J. Saxe, Embeddability of weighted graphs in $k$-space is strongly NP-hard, in: Proceedings of 17th Allerton Conference in Communications, Control and Computing, 1979, pp. 480–489.

[43] T. Schlick, Molecular Modelling and Simulation: An Interdisciplinary Guide, Springer, New York, 2002.

[44] I. Schoenberg, Remarks to Maurice Fréchet's article Sur la définition axiomatique d'une classe d'espaces distanciés vectoriellement applicable sur l'espace de Hilbert, Annals of Mathematics 36 (3) (1935) 724–732.

[45] A.M.-C. So, Y. Ye, Theory of semidefinite programming for sensor network localization, Mathematical Programming, Series B 109 (2007) 367–384.

[46] T.-S. Tay, W. Whiteley, Generating isostatic frameworks, Structural Topology 11 (1985) 21–69.

[47] W. Whiteley, Infinitesimally rigid polyhedra. I. Statics of frameworks, Transactions of the American Mathematical Society 285 (2) (1984) 431–465.

[48] W. Whiteley, Infinitesimally rigid polyhedra. II: modified spherical frameworks, Transactions of the American Mathematical Society 306 (1) (1988) 115–139.

[49] Z. Zhu, A.M.-C. So, Y. Ye, Universal rigidity and edge sparsification for sensor network localization, SIAM Journal on Optimization 20 (6) (2010) 3059–3081.