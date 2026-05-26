# Scopus Screening Batch 01

Records: 20
Result file to write: `references/screening/agent-results/scopus-batch-01-screening.md`

Classify each record as `core`, `background`, `maybe`, or `reject`.
Use the tag vocabulary from the reference-screen skill. Give one sentence
of rationale for every `core` and `maybe` record; terse rationale is enough
for `background` and `reject` records.

## R001. An impossible combinatorial counting method in distance geometry

- Year: 2024
- Venue: Discrete Applied Mathematics
- DOI: 10.1016/j.dam.2024.02.018
- Cited by: 1
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts, s4-ddgp-counting-limits, s5-lateration
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, solution counting, distance geometry applications
- Authors: Abud G.; Alencar J.; Lavor C.; Liberti L.; Mucherino A.
- Author keywords: Branch-and-Prune; DDGP; DGP; DMDGP; Solution symmetry
- Index keywords: Graph theory; Branch and prunes; DDGP; DGP; Distance geometry; DMDGP; Fixed-parameter tractability; Geometric representation; Geometry problems; Solution symmetry; Vertex ordering; Geometry

Abstract: The Distance Geometry Problem asks for a geometric representation of a given weighted graph in RK so that vertices are points and edges are segments with lengths equal to the corresponding weights. Two problem variants are defined by a vertex order given as part of the input, which allows the use of a branching algorithm based on K-lateration: find two possible positions for the next vertex j using the positions of K predecessors and their distances to j, then explore each position recursively, pruning positions at need. Whereas the first variant only requires the K predecessors to exist, the second variant also requires them to be contiguous and immediately preceding j. For this variant, fixed-parameter tractability of the algorithm can be established by means of a solution counting method that only depends on the graph edges (rather than their weights). Only in the first variant, however, it is possible to efficiently construct a suitable vertex order directly from the graph. Since both fixed-parameter tractability and efficient vertex order construction are desirable properties, ...

## R002. On the optimality of finding DMDGP symmetries

- Year: 2021
- Venue: Computational and Applied Mathematics
- DOI: 10.1007/s40314-021-01479-6
- Cited by: 4
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, distance geometry applications
- Authors: Lavor C.; Oliveira A.; Rocha W.; Souza M.
- Author keywords: Branch-and-prune; DMDGP; DMDGP symmetries; Protein geometry
- Index keywords: Geometry; Branch and prunes; Discretizable molecular distance geometry problem; Discretizable molecular distance geometry problem symmetry; Distance geometry; Geometry problems; Molecular distance geometry problem; Optimality; Protein geometry; Simple++; Undirected graphs

Abstract: The Discretizable Molecular Distance Geometry Problem (DMDGP) is a subclass of the Distance Geometry Problem, which aims to embed a weighted simple undirected graph in a Euclidean space, such that the distances between the points correspond to the values given by the weighted edges in the graph. The search space of the DMDGP is combinatorial, based on a total vertex order that implies symmetry properties related to partial reflections around planes defined by the Cartesian coordinates of three immediate and consecutive vertices that precede the so-called symmetry vertices. Since these symmetries allow us to know a priori the cardinality of the solution set and to calculate all the DMDGP solutions, given one of them, it would be relevant to identify these symmetries efficiently. Exploiting mathematical properties of the vertices associated with these symmetries, we present an optimal algorithm that finds such vertices. © 2021, SBMAC - Sociedade Brasileira de Matemática Aplicada e Computacional.

## R003. Exploiting symmetry properties of the discretizable molecular distance geometry problem

- Year: 2012
- Venue: Journal of Bioinformatics and Computational Biology
- DOI: 10.1142/S0219720012420097
- Cited by: 27
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, solution counting, distance geometry applications
- Authors: Mucherino A.; Lavor C.; Liberti L.
- Author keywords: branch-and-prune; discretization; distance geometry; NMR; Protein conformation; symmetry
- Index keywords: Algorithms; Models, Molecular; Proteins; protein; algorithm; article; chemical structure; chemistry

Abstract: The Discretizable Molecular Distance Geometry Problem (DMDGP) involves a subset of instances of the distance geometry problem for which some assumptions allowing for discretization are satisfied. The search domain for the DMDGP is a binary tree that can be efficiently explored by employing a Branch & Prune (BP) algorithm. We showed in recent works that this binary tree may contain several symmetries, which are directly related to the total number of solutions of DMDGP instances. In this paper, we study the possibility of exploiting these symmetries for speeding up the solution of DMDGPs, and propose an extension of the BP algorithm that we named symmetry-driven BP (symBP). Computational experiments on artificial and protein instances are presented. © 2012 Imperial College Press.

## R004. On the number of solutions of the discretizable molecular distance geometry problem

- Year: 2011
- Venue: Lecture Notes in Computer Science (including subseries Lecture Notes in Artificial Intelligence and Lecture Notes in Bioinformatics)
- DOI: 10.1007/978-3-642-22616-8_26
- Cited by: 12
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, solution counting, distance geometry applications
- Authors: Liberti L.; Masson B.; Lee J.; Lavor C.; Mucherino A.
- Author keywords: Branch-and-Prune; distance geometry; power of two; symmetry
- Index keywords: Combinatorial optimization; Crystal symmetry; Branch and prunes; Combinatorial algorithm; Distance geometry; Molecular distance geometry problem; Power-of-two; Geometry

Abstract: The Discretizable Molecular Distance Geometry Problem is a subset of instances of the distance geometry problem that can be solved by a combinatorial algorithm called "Branch-and-Prune". It was observed empirically that the number of solutions of YES instances is always a power of two. We perform an extensive theoretical analysis of the number of solutions for these instances and we prove that this number is a power of two with probability one. © 2011 Springer-Verlag.

## R005. The discretizable distance geometry problem

- Year: 2012
- Venue: Optimization Letters
- DOI: 10.1007/s11590-011-0358-3
- Cited by: 64
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Mucherino A.; Lavor C.; Liberti L.
- Author keywords: Branch and prune; Combinatorial reformulations; DDGP<sub>3</sub>; Distance geometry; DMDGP
- Index keywords: Proteins; Branch and prune; Combinatorial reformulations; ; Distance geometry; DMDGP; Geometry

Abstract: We introduce the discretizable distance geometry problem in ℝ3 (DDGP3), which consists in a subclass of instances of the Distance Geometry Problem for which an embedding in ℝ3 can be found by means of a discrete search. We show that the DDGP3 is a generalization of the discretizable molecular distance geometry problem (DMDGP), and we discuss the main differences between the two problems. We prove that the DDGP3 is NP-hard and we extend the Branch & Prune (BP) algorithm, previously used for the DMDGP, for solving instances of the DDGP3. Protein graphs may or may not be in DMDGP and/or DDGP3 depending on vertex orders and edge density. We show experimentally that as distance thresholds decrease, PDB protein graphs which fail to be in the DMDGP still belong to DDGP3, which means that they can still be solved using a discrete search. © 2011 Springer-Verlag.

## R006. A New Algorithm for the K DMDGP Subclass of Distance Geometry Problems with Exact Distances

- Year: 2021
- Venue: Algorithmica
- DOI: 10.1007/s00453-021-00835-6
- Cited by: 4
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, distance geometry applications
- Authors: Gonçalves D.S.; Lavor C.; Liberti L.; Souza M.
- Author keywords: Discretization; Distance Geometry; DMDGP; Symmetries
- Index keywords: Binary trees; Geometry; Trees (mathematics); Branch and prunes; Combinatorial search; Computational results; Distance geometry; Molecular distance geometry problem; Point distances; Protein conformation; Symmetry properties; Inverse problems

Abstract: The fundamental inverse problem in distance geometry is the one of finding positions from inter-point distances. The Discretizable Molecular Distance Geometry Problem (DMDGP) is a subclass of the Distance Geometry Problem (DGP) whose search space can be discretized and represented by a binary tree, which can be explored by a Branch-and-Prune (BP) algorithm. It turns out that this combinatorial search space possesses many interesting symmetry properties that were studied in the last decade. In this paper, we present a new algorithm for this subclass of the DGP, which exploits DMDGP symmetries more effectively than its predecessors. Computational results show that the speedup, with respect to the classic BP algorithm, is considerable for sparse DMDGP instances related to protein conformation. © 2021, The Author(s), under exclusive licence to Springer Science+Business Media, LLC, part of Springer Nature.

## R007. The discretizable molecular distance geometry problem seems easier on proteins

- Year: 2013
- Venue: Distance Geometry: Theory, Methods, and Applications
- DOI: 10.1007/978-1-4614-5128-0_3
- Cited by: 22
- Document type: Book chapter
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, distance geometry applications
- Authors: Liberti L.; Lavor C.; Mucherino A.
- Author keywords: Branch-and-Prune; Distance geometry; Fixed-parameter tractable; Protein conformation; Symmetry
- Index keywords: Geometry; Integer programming; Nuclear magnetic resonance; Branch and prunes; Branch-and-prune algorithms; Distance geometry; Embeddings; Fixed-parameter tractable; Inter-atomic distances; Molecular conformation; Molecular distance geometry problem; Protein conformation; Symmetry; Proteins

Abstract: Distance geometry methods are used to turn a set of interatomic distances given by Nuclear Magnetic Resonance (NMR) experiments into a consistent molecular conformation. In a set of papers (see the survey [8]) we proposed a Branch-and-Prune (BP) algorithm for computing the set X of all incongruent embeddings of a given protein backbone. Although BP has a worst-case exponential running time in general, we always noticed a linear-like behaviour in computational experiments. In this chapter we provide a theoretical explanation to our observations. We show that the BP is fixed-parameter tractable on protein-like graphs and empirically show that the parameter is constant on a set of proteins from the Protein Data Bank. © 2013 Springer Science+Business Media New York. All rights reserved.

## R008. Geometric Algebra to Describe the Exact Discretizable Molecular Distance Geometry Problem for an Arbitrary Dimension

- Year: 2019
- Venue: Advances in Applied Clifford Algebras
- DOI: 10.1007/s00006-019-0995-7
- Cited by: 8
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, distance geometry applications, algebraic/genericity methods
- Authors: Camargo V.S.; Castelani E.V.; Fernandes L.A.F.; Fidalgo F.
- Author keywords: Branch-and-prune; Conformal algebra; Distance geometry problem; Minkowski space; Outer product; Symmetry
- Index keywords: n/a

Abstract: The K-Discretizable Molecular Distance Geometry Problem (KDMDGP) is a subclass of the Distance Geometry Problem (DGP), whose complexity is NP-hard, such that the search space is finite. In this work, the authors describe it completely using Conformal Geometric Algebra (CGA), exploring a Minkowski space that provides a natural interpretation of hyperspheres, hyperplanes, points and pair of points as computational primitives, which are largely relevant to the KDMDGP. It also presents a theoretical approach to solve the KDMDGP using ideas from classic Branch-and-Prune (BP) algorithm in this new fashion. Time complexity analysis and practical computational results showed that the naive implementation of the CGA is not as efficient as classical formulation. In order to illustrate this, preliminary results are displayed at the end and, also, directions to future developments. © 2019, Springer Nature Switzerland AG.

## R009. A Probabilistic Approach in the Search Space of the Molecular Distance Geometry Problem

- Year: 2025
- Venue: Journal of Chemical Information and Modeling
- DOI: 10.1021/acs.jcim.4c00427
- Cited by: 0
- Document type: Article
- Query hits: s1-core-ddgp, s3-symmetry-counts, s6-local-symmetry
- Tags hint: DDGP/DMDGP definitions, partial reflection symmetries, distance geometry applications
- Authors: Marques R.S.; Souza M.; Batista F.; Gonçalves M.; Lavor C.
- Author keywords: n/a
- Index keywords: Algorithms; Databases, Protein; Magnetic Resonance Spectroscopy; Models, Molecular; Nuclear Magnetic Resonance, Biomolecular; Probability; Protein Conformation; Proteins; Geometry; Molecular docking; Nuclear magnetic resonance; Trees (mathematics); protein; Binary string; Depth first; Distance information; Inter-atomic distances; Molecular distance geometry problem; Probabilistics approach; Protein data bank; Protein molecules; Search spaces; Three-dimensional shape; algorithm; chemistry; heteronuclear nuclear magnetic resonance; molecular model; nuclear magnetic resonance spectroscopy; probability; procedures; protein conformation; protein database; Binary trees

Abstract: The discovery of the three-dimensional shape of protein molecules using interatomic distance information from nuclear magnetic resonance (NMR) can be modeled as a discretizable molecular distance geometry problem (DMDGP). Due to its combinatorial characteristics, the problem is conventionally solved in the literature as a depth-first search in a binary tree. In this work, we introduce a new search strategy, which we call frequency-based search (FBS), that for the first time utilizes geometric information contained in the protein data bank (PDB). We encode the geometric configurations of 14,382 molecules derived from NMR experiments present in the PDB into binary strings. The obtained results show that the sample space of the binary strings extracted from the PDB does not follow a uniform distribution. Furthermore, we compare the runtime of the symmetry-based build-Up (SBBU) algorithm (the most efficient method in the literature to solve the DMDGP) combined with FBS and the depth-first search (DFS) in finding a solution, ascertaining that FBS performs better in about 70% of the ...

## R010. A symmetry-driven BP algorithm for the discretizable molecular distance geometry problem

- Year: 2011
- Venue: 2011 IEEE International Conference on Bioinformatics and Biomedicine Workshops, BIBMW 2011
- DOI: 10.1109/BIBMW.2011.6112403
- Cited by: 5
- Document type: Conference paper
- Query hits: s1-core-ddgp, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, partial reflection symmetries, solution counting, distance geometry applications
- Authors: Mucherino A.; Lavor C.; Liberti L.
- Author keywords: n/a
- Index keywords: Binary trees; BP algorithm; Computational experiment; Deterministic algorithms; Molecular distance geometry problem; NMR techniques; Structural biology; Backpropagation algorithms

Abstract: Branch & Prune (BP) is a deterministic algorithm for the solution of the Discretizable Molecular Distance Geometry Problem (DMDGP). This problem has important applications in the field of structural biology, in particular for the determination of the three-dimensional conformation of a molecule by using information obtained by NMR techniques. In recent works, we proved that the search domain of the DMDGP, which is represented by a binary tree, contains various symmetries which are related to the number of solutions to the problem. In the present work, we propose a variant of the BP algorithm which is able to exploit the information regarding the symmetries to speed up the search. Computational experiments show that the symmetry-driven BP (symBP) outperforms the original BP algorithm in particular when instances having several solutions are considered. © 2011 IEEE.

## R011. The discretizable molecular distance geometry problem

- Year: 2012
- Venue: Computational Optimization and Applications
- DOI: 10.1007/s10589-011-9402-6
- Cited by: 98
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Lavor C.; Liberti L.; MacUlan N.; Mucherino A.
- Author keywords: Branch-and-prune; Distance geometry; Molecular conformation; NMR; Protein
- Index keywords: Computational methods; Nuclear magnetic resonance; Optimization; Proteins; BP algorithm; Branch-and-prune; Computational results; Discrete spaces; Distance geometry; Molecular conformation; Molecular distance geometry problem; NP-hard; Solution accuracy; Solution algorithms; Weighted undirected graph; Geometry

Abstract: Given a simple weighted undirected graph G=(V,E,d) with d:E→ℝ +, the Molecular Distance Geometry Problem (MDGP) consists in finding an embedding x:V→ℝ 3 such that ∥x u -x v ∥=d uv for each {u,v} E. We show that under a few assumptions usually satisfied in proteins, the MDGP can be formulated as a search in a discrete space. We call this MDGP subclass the Discretizable MDGP (DMDGP). We show that the DMDGP is NP-hard and we propose a solution algorithm called Branch-and-Prune (BP). The BP algorithm performs remarkably well in practice in terms of speed and solution accuracy, and can be easily modified to find all incongruent solutions to a given DMDGP instance. We show computational results on several artificial and real-life instances. © 2011 Springer Science+Business Media, LLC.

## R012. Discretization vertex orders in distance geometry

- Year: 2015
- Venue: Discrete Applied Mathematics
- DOI: 10.1016/j.dam.2014.08.035
- Cited by: 42
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Cassioli A.; Günlük O.; Lavor C.; Liberti L.
- Author keywords: Branch-and-Prune; DDGP; DMDGP; Molecular conformation; Protein; Re-order
- Index keywords: Mathematical programming; Proteins; Branch and prunes; DDGP; DMDGP; Molecular conformation; Re orderings; Graph algorithms

Abstract: When a weighted graph is an instance of the Distance Geometry Problem (DGP), certain types of vertex orders (called discretization orders) allow the use of a very efficient, precise and robust discrete search algorithm (called Branch-and-Prune). Accordingly, finding such orders is critically important in order to solve DGPs in practice. We discuss three types of discretization orders, the complexity of determining their existence in a given graph, and the inclusion relations between the three order existence problems. We also give three mathematical programming formulations of some of these ordering problems. © 2014 Elsevier B.V. All rights reserved.

## R013. A symmetry-based splitting strategy for discretizable distance geometry problems

- Year: 2018
- Venue: Journal of Global Optimization
- DOI: 10.1007/s10898-018-0610-9
- Cited by: 8
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, distance geometry applications
- Authors: Fidalgo F.; Gonçalves D.S.; Lavor C.; Liberti L.; Mucherino A.
- Author keywords: Branch-and-prune; Decomposition; Partial reflection; Protein structure determination
- Index keywords: Decomposition; Forestry; Trees (mathematics); Algorithm performance; Branch and prunes; Branch-and-prune algorithms; Computational experiment; Exhaustive enumeration; Partial reflection; Protein structures; Splitting strategies; Geometry

Abstract: Discretizable distance geometry problems consist in a subclass of distance geometry problems where the search space can be discretized and reduced to a tree. Such problems can be tackled by applying a branch-and-prune algorithm, which is able to perform an exhaustive enumeration of the solution set. In this work, we exploit the concept of symmetry in the search tree for isolating subtrees that are explored only one time for improving the algorithm performances. The proposed strategy is based on the idea of dividing an original instance of the problem into sub-instances that can thereafter be solved (almost) independently. We present some computational experiments on a set of artificially generated instances, with exact distances, to validate the theoretical results. © 2018, Springer Science+Business Media, LLC, part of Springer Nature.

## R014. Geometric Algebra to Model Uncertainties in the Discretizable Molecular Distance Geometry Problem

- Year: 2017
- Venue: Advances in Applied Clifford Algebras
- DOI: 10.1007/s00006-016-0653-2
- Cited by: 27
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications, algebraic/genericity methods
- Authors: Alves R.; Lavor C.
- Author keywords: 3D protein structure; Branch and prune algorithm; Conformal geometric algebra; Distance geometry
- Index keywords: n/a

Abstract: The discretizable molecular distance geometry problem (DMDGP) is related to the determination of 3D protein structure using distance information detected by nuclear magnetic resonance (NMR) experiments. The chemistry of proteins and the NMR distance information allow us to define an atomic order v1, … , vn such that the distances related to the pairs { vi - 3, vi} , { vi - 2, vi} , { vi - 1, vi} , for i> 3 , are available, which implies that the search space can be represented by a tree. A DMDGP solution can be represented by a path from the root to a leaf node of this tree, found by an exact method, called branch-and-prune (BP). Because of uncertainty in NMR data, some of the distances related to the pairs { vi - 3, vi} may not be precise values, being represented by intervals of real numbers [d̲i-3,i,d¯i-3,i]. In order to apply BP algorithm in this context, sample values from those intervals should be taken. The main problem of this approach is that if we sample many values, the search space increases drastically, and for small samples, no solution can be found. We explain how ...

## R015. Tuning interval Branch-and-Prune for protein structure determination

- Year: 2018
- Venue: Journal of Global Optimization
- DOI: 10.1007/s10898-018-0635-0
- Cited by: 20
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Worley B.; Delhommel F.; Cordier F.; Malliavin T.E.; Bardiaux B.; Wolff N.; Nilges M.; Lavor C.; Liberti L.
- Author keywords: Branch-and-Prune; Distance geometry; Nuclear magnetic resonance; Protein structure
- Index keywords: Application programs; Iterative methods; Nuclear magnetic resonance; Proteins; Branch and prunes; Distance geometry; Existing structure; Molecular distance geometry problem; Orders of magnitude; Protein structures; Research and application; Structure determination; Bioinformatics

Abstract: The interval Branch and Prune (iBP) algorithm for obtaining solutions to the interval Discretizable Molecular Distance Geometry Problem (iDMDGP) has proven itself as a powerful method for molecular structure determination. However, substantial obstacles still must be overcome before iBP may be employed as a tractable general-purpose alternative to existing structure determination algorithms. This work introduces an iterative variant of the iBP algorithm that leverages existing knowledge of protein structures in order to reduce the size of the effective search space by many orders of magnitude. These improvements are included in a newly released implementation of the iBP software that aims to provide a solid platform for both research and application of the iDMDGP. © 2018, Springer Science+Business Media, LLC, part of Springer Nature.

## R016. The interval Branch-and-Prune algorithm for the discretizable molecular distance geometry problem with inexact distances

- Year: 2013
- Venue: Journal of Global Optimization
- DOI: 10.1007/s10898-011-9799-6
- Cited by: 77
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Lavor C.; Liberti L.; Mucherino A.
- Author keywords: Combinatorial optimization; Distance geometry; Interval Branch and Prune; NMR data; Protein conformations
- Index keywords: Algorithms; Combinatorial optimization; Geometry; Proteins; Branch and prunes; Branch-and-prune algorithms; Computational experiment; Distance geometry; Euclidean distance; Molecular distance geometry problem; NMR data; Protein conformation; Problem solving

Abstract: The Distance Geometry Problem in three dimensions consists in finding an embedding in {\mathbb R3} of a given nonnegatively weighted simple undirected graph such that edge weights are equal to the corresponding Euclidean distances in the embedding. This is a continuous search problem that can be discretized under some assumptions on the minimum degree of the vertices. In this paper we discuss the case where we consider the full-atom representation of the protein backbone and some of the edge weights are subject to uncertainty within a given nonnegative interval. We show that a discretization is still possible and propose the iBP algorithm to solve the problem. The approach is validated by some computational experiments on a set of artificially generated instances. © 2011 Springer Science+Business Media, LLC.

## R017. Recent advances on the Discretizable Molecular Distance Geometry Problem

- Year: 2012
- Venue: European Journal of Operational Research
- DOI: 10.1016/j.ejor.2011.11.007
- Cited by: 61
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Lavor C.; Liberti L.; MacUlan N.; Mucherino A.
- Author keywords: Bioinformatics; Branch-and-prune; Graph theory; Protein conformation
- Index keywords: Bioinformatics; Graph theory; Molecular biology; Adjacent vertices; Branch-and-prune; Continuous spaces; Edge weights; Euclidean distance; Molecular distance geometry problem; NMR data; Protein conformation; Protein structures; Proteomics; Undirected graph; Proteins

Abstract: The Molecular Distance Geometry Problem (MDGP) consists in finding an embedding in R 3 of a nonnegatively weighted simple undirected graph with the property that the Euclidean distances between embedded adjacent vertices must be the same as the corresponding edge weights. The Discretizable Molecular Distance Geometry Problem (DMDGP) is a particular subset of the MDGP which can be solved using a discrete search occurring in continuous space; its main application is to find three-dimensional arrangements of proteins using Nuclear Magnetic Resonance (NMR) data. The model provided by the DMDGP, however, is too abstract to be directly applicable in proteomics. In the last five years our efforts have been directed towards adapting the DMDGP to be an ever closer model of the actual difficulties posed by the problem of determining protein structures from NMR data. This survey lists recent developments on DMDGP related research. © 2011 Elsevier B.V. All rights reserved.

## R018. Branch-and-prune trees with bounded width

- Year: 2011
- Venue: 10th Cologne-Twente Workshop on Graphs and Combinatorial Optimization, CTW 2011 - Proceedings of the Conference
- DOI: n/a
- Cited by: 7
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, partial reflection symmetries, distance geometry applications
- Authors: Liberti L.; Masson B.; Lavor C.; Mucherino A.
- Author keywords: Distance geometry; DMDGP; Order; Reflection; Symmetry
- Index keywords: n/a

Abstract: [No abstract available]

## R019. Solving the discretizable molecular distance geometry problem by multiple realization trees

- Year: 2013
- Venue: Distance Geometry: Theory, Methods, and Applications
- DOI: 10.1007/978-1-4614-5128-0_9
- Cited by: 5
- Document type: Book chapter
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Nucci P.; Nogueira L.T.; Lavor C.
- Author keywords: Branch-and-prune; Distance geometry; Realization tree
- Index keywords: Geometry; Trees (mathematics); Branch and prunes; Branch-and-prune algorithms; Discrete algorithms; Distance geometry; Molecular distance geometry problem; Realization tree; Search spaces; Forestry

Abstract: The discretizable molecular distance geometry problem (DMDGP) is a subclass of the MDGP, where instances can be solved using a discrete algorithm called branch-and-prune (BP). We present an initial study showing that the BP algorithm can be used differently from its original form, where the initial atoms are fixed and the branches of the BP tree are generated until the last atom is reached. Particularly, we show that the use of multiple BP trees may explore the search space faster than the original BP. © 2013 Springer Science+Business Media New York. All rights reserved.

## R020. Advances on the geometric algebra approach to the discretizable molecular distance geometry problem (DMDGP)

- Year: 2016
- Venue: ACM International Conference Proceeding Series
- DOI: 10.1145/2949035.2949057
- Cited by: 2
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications, algebraic/genericity methods
- Authors: Alves R.; De Souza C.; Lavor C.
- Author keywords: 3D protein structure; Branch and prune algorithm; Conformal geometric algebra; Distance geometry
- Index keywords: Algebra; Bioinformatics; Computer graphics; Nuclear magnetic resonance; Proteins; Atomic distances; Branch-and-prune algorithms; Conformal Geometric Algebra; Distance geometry; Geometric Algebra; Molecular distance geometry problem; Nuclear magnetic resonance(NMR); Protein structures; Geometry

Abstract: We consider the Discretizable Molecular Distance Geometry Problem (DMDGP), which is a subclab of the Distance Geometry Problem (DGP), where the search space can be discretized. The DMDGP consists in finding a 3D protein structure for which some of the atomic distances are provided by Nuclear Magnetic Resonance (NMR) experiments. Due to the uncertainties in NMR data, some distances may not be precise, being represented by intervals of real numbers. We present recent results related to the use of Conformal Geometric Algebra to model the DMDGP with imprecise data. © 2016 ACM.
