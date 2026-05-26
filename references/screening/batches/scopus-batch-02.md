# Scopus Screening Batch 02

Records: 20
Result file to write: `references/screening/agent-results/scopus-batch-02-screening.md`

Classify each record as `core`, `background`, `maybe`, or `reject`.
Use the tag vocabulary from the reference-screen skill. Give one sentence
of rationale for every `core` and `maybe` record; terse rationale is enough
for `background` and `reject` records.

## R021. Constraint programming approaches for the discretizable molecular distance geometry problem

- Year: 2022
- Venue: Networks
- DOI: 10.1002/net.22068
- Cited by: 0
- Document type: Article
- Query hits: s1-core-ddgp, s3-symmetry-counts, s5-lateration
- Tags hint: DDGP/DMDGP definitions, partial reflection symmetries, distance geometry applications
- Authors: MacNeil M.; Bodur M.
- Author keywords: combinatorial optimization; constraint programming; contiguous trilateration ordering; discretizable molecular distance geometry; discretization order
- Index keywords: Computer programming; Constraint theory; Geometry; Graph structures; Computational results; Constraint programming; Distance geometry; Integer programming formulations; Molecular distance geometry problem; Reduction techniques; State of the art; Symmetry breaking constraints; Integer programming

Abstract: The Distance Geometry Problem (DGP) seeks to find positions for a set of points in geometric space when some distances between pairs of these points are known. The so-called discretization assumptions allow us to discretize the search space of DGP instances. In this paper, we focus on a key subclass of DGP, namely the Discretizable Molecular DGP, and study its associated graph vertex ordering problem, the Contiguous Trilateration Ordering Problem (CTOP), which helps solve DGP. We propose the first constraint programming formulations for CTOP, as well as a set of checks for proving infeasibility, domain reduction techniques, symmetry breaking constraints, and valid inequalities. Our computational results on random and pseudo-protein instances indicate that our formulations outperform the state-of-the-art integer programming formulations. © 2021 Wiley Periodicals LLC.

## R022. Comparisons between an exact and a metaheuristic algorithm for the molecular distance geometry problem

- Year: 2009
- Venue: Proceedings of the 11th Annual Genetic and Evolutionary Computation Conference, GECCO-2009
- DOI: 10.1145/1569901.1569948
- Cited by: 23
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Mucherino A.; Liberti L.; Lavor C.; Maculan N.
- Author keywords: Branch and prune; Combinatorial optimization; Distance geometry; Monkey search; Protein molecules
- Index keywords: Combinatorial optimization; Conformations; Food supply; Genetic algorithms; Geometry; Nuclear magnetic resonance; Proteins; Trees (mathematics); Branch and prunes; Distance geometry; Meta heuristic algorithm; Molecular conformation; Molecular distance geometry problem; Monkey search; Nuclear Magnetic Resonance (NMR); Protein molecules; Heuristic algorithms

Abstract: We consider the Discretizable Molecular Distance Geometry Problem (DMDGP), which consists in a subclass of instances of the distance geometry problem related to molecular conformations for which a combinatorial reformulation can be supplied. We investigate the performances of two different algorithms for solving the DMDGP. The first one is the Branch and Prune (BP) algorithm, an exact algorithm that is strongly based on the structure of the combinatorial problem. The second one is the Monkey Search (MS) algorithm, a meta-heuristic algorithm that is inspired by the behavior of a monkey climbing trees in search for food supplies, and that exploits ideas and strategies from other meta-heuristic searches, such Genetic Algorithms, Differential Evolution, and so on. The comparison between the two algorithms is performed on a set of instances related to protein conformations. The used instances simulate data obtained from the Nuclear Magnetic Resonance (NMR), because the typical distances provided by NMR are considered and a predetermined number of wrong distances are included. Copyright ...

## R023. A parallel BP algorithm for the discretizable distance geometry problem

- Year: 2012
- Venue: Proceedings of the 2012 IEEE 26th International Parallel and Distributed Processing Symposium Workshops, IPDPSW 2012
- DOI: 10.1109/IPDPSW.2012.218
- Cited by: 5
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Gramacho W.; Mucherino A.; Lavor C.; Maculan N.
- Author keywords: branch-and-prune; discretization; distance geometry; parallel computing
- Index keywords: Distributed parameter networks; Geometry; Parallel algorithms; Parallel architectures; BP algorithm; branch-and-prune; Common coordinate systems; Computational experiment; Discretizations; Distance geometry; Parallel version; Vertex set; Parallel processing systems

Abstract: We propose a parallel version of the Branch & Prune (BP) algorithm for the Discretizable Distance Geometry Problem (DDGP), which consists in a subclass of Distance Geometry Problems (DGPs) that can be discretized. The main idea is to split a DDGP instance in as many sub instances as the number of processors involved in the computation, and to invoke the sequential version of BP on each processor. Due to the flexibility of the discretizing orderings that can be defined on the vertex sets of graphs G representing DDGP instances, the subdivision of the original instance can be performed so that all solutions generated by locally solving the several sub instances are represented in a common coordinate system. This way, the communication phase of the parallel algorithm, where the local solutions are combined in order to generate the final set of solutions, is very efficient. We present some preliminary computational experiments and we study the behavior of the algorithm in relation to the number of considered processors. We also give some directions for transforming DDGP instances in ...

## R024. An efficient exhaustive search for the discretizable distance geometry problem with interval data

- Year: 2019
- Venue: Proceedings of the 2019 Federated Conference on Computer Science and Information Systems, FedCSIS 2019
- DOI: 10.15439/2019F62
- Cited by: 5
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Mucherino A.; Lin J.-H.
- Author keywords: n/a
- Index keywords: Binary trees; Forestry; Information systems; Information use; Parameter estimation; Computational experiment; Distance constraints; Distance geometry; Distance information; Initial estimation; Local optimizations; Resolution parameters; Weighted undirected graph; Trees (mathematics)

Abstract: The Distance Geometry Problem (DGP) asks whether a simple weighted undirected graph can be realized in a given space (generally Euclidean) so that a given set of distance constraints (associated to the edges of the graph) is satisfied. The Discretizable DGP (DDGP) represents a subclass of instances where the search space can be reduced to a discrete domain having the structure of a tree. In the ideal case where all distances are precise, the tree is binary and one singleton, representing one possible position for a vertex of the graph, is associated to every tree node. When the distance information is however not precise, the uncertainty on the distance values implies that a three-dimensional region of the search space needs to be assigned to some nodes of the tree. By using a recently proposed coarse-grained representation for DDGP solutions, we extend in this work the branch-and-prune (BP) algorithm so that it can efficiently perform an exhaustive search of the search domain, even when the uncertainty on the distances is important. Instead of associating singletons to nodes, we ...

## R025. On the number of realizations of certain Henneberg graphs arising in protein conformation

- Year: 2014
- Venue: Discrete Applied Mathematics
- DOI: 10.1016/j.dam.2013.01.020
- Cited by: 45
- Document type: Article
- Query hits: s2-branch-prune
- Tags hint: Branch-and-Prune, partial reflection symmetries, solution counting, distance geometry applications
- Authors: Liberti L.; Masson B.; Lee J.; Lavor C.; Mucherino A.
- Author keywords: Branch-and-Prune; Distance geometry; Graph rigidity; Partial reflection; Protein conformation
- Index keywords: Rigidity; Adjacent vertices; Application fields; Branch and prunes; Distance geometry; Edge-weighted graph; Euclidean coordinates; Partial reflection; Protein conformation; Proteins

Abstract: Several application fields require finding Euclidean coordinates consistent with a set of distances. More precisely, given a simple undirected edge-weighted graph, we wish to find a realization in a Euclidean space so that adjacent vertices are placed at a distance which is equal to the corresponding edge weight. Realizations of a graph can be either flexible or rigid. In certain cases, rigidity can be seen as a property of the graph rather than the realization. In the last decade, several advances have been made in graph rigidity, but most of these, for applicative reasons, focus on graphs having a unique realization. In this paper we consider a particular type of weighted Henneberg graphs that model protein backbones and show that almost all of them give rise to sets of incongruent realizations whose cardinality is a power of two. © 2013 Elsevier B.V. All rights reserved.

## R026. An autonomic parallel strategy for exhaustive search tree algorithms on shared or heterogeneous systems

- Year: 2024
- Venue: Concurrency and Computation: Practice and Experience
- DOI: 10.1002/cpe.7955
- Cited by: 5
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Passos F.G.O.; Rebello V.E.F.
- Author keywords: autonomic computing; heterogeneous distributed systems; parallel branch-and-prune algorithms; search tree algorithms; self-configuring parallel algorithms
- Index keywords: Computing power; Trees (mathematics); Autonomic Computing; Autonomics; Branch-and-prune algorithms; Heterogeneous distributed systems; Molecular distance geometry problem; Parallel branch-and-prune algorithm; Parallel branches; Search tree algorithms; Self-configuring; Self-configuring parallel algorithm; Middleware

Abstract: Backtracking branch-and-prune (BP) algorithms and their variants are exhaustive search tree techniques widely employed to solve optimization problems in many scientific areas. However, they characteristically often demand significant amounts of computing power for problem sizes representative of real-world scenarios. Given that their search domains can often be partitioned, these algorithms are frequently designed to execute in parallel by harnessing distributed computing systems. However, to achieve efficient parallel execution times, an effective strategy is required to balance the nonuniform partition workloads across the available resources. Furthermore, with the increasing integration of servers with heterogeneous resources and the adoption of resource sharing, balancing workloads is becoming complex. This paper proposes a strategy to execute parallel BP algorithms more efficiently on even shared or heterogeneous distributed systems. The approach integrates a self-adjusting dynamic partitioning method in the BP algorithm with a dynamic scheduler, provided by an application ...

## R027. A numerical-and-computational study on the impact of using quaternions in the branch-and-prune algorithm for exact discretizable distance geometry problems

- Year: 2024
- Venue: Computational Optimization and Applications
- DOI: 10.1007/s10589-023-00526-8
- Cited by: 0
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Fidalgo F.; Castelani E.; Philippi G.
- Author keywords: Branch-and-prune; Distance geometry; Numerical optimization; Protein conformation; Quaternions
- Index keywords: Computational geometry; Global optimization; Integer programming; Proteins; Branch and prunes; Branch-and-prune algorithms; Computational studies; Distance geometry; Geometric relations; Geometry problems; Numerical optimizations; Primitive element; Protein conformation; Quaternion; Binary trees

Abstract: Distance geometry is a branch of Mathematics which studies geometric relations having distances as primitive elements. Its fundamental problem, the distance geometry problem, consists in determining point positions in RK such that their Euclidean distances match those in a given list of inter-point distances. Such problem can be cast as a global optimization problem which is often tackled with continuous optimization techniques. If K=3, it is called molecular DGP (MDGP). Under assumptions on the available distances in this list, the search space for the MDGP can be discretized so that it is able to be designed as a binary tree, giving birth to the discretized MDGP. The main method to solve it is the Branch-and-Prune Algorithm, a recursive and combinatorial tool that explores such tree in a depth-first search and whose classical formulation is based in a homogeneous matrix product that encodes one translation and two rotations. This paper presents a numerical study on the theoretical computational effort to do that task together with a quaternionic proposal as an alternative for the ...

## R028. A Coarse-Grained Representation for Discretizable Distance Geometry with Interval Data

- Year: 2019
- Venue: Lecture Notes in Computer Science (including subseries Lecture Notes in Artificial Intelligence and Lecture Notes in Bioinformatics)
- DOI: 10.1007/978-3-030-17938-0_1
- Cited by: 4
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Mucherino A.; Lin J.-H.; Gonçalves D.S.
- Author keywords: n/a
- Index keywords: Bioinformatics; Biomedical engineering; Forestry; Branch and prunes; Computational experiment; Distance constraints; Distance geometry; Distance information; Inter-atomic distances; Projected gradient; Real-life applications; Atoms

Abstract: We propose a coarse-grained representation for the solutions of discretizable instances of the Distance Geometry Problem (DGP). In several real-life applications, the distance information is not provided with high precision, but an approximation is rather given. We focus our attention on protein instances where inter-atomic distances can be either obtained from the chemical structure of the molecule (which are exact), or through experiments of Nuclear Magnetic Resonance (which are generally represented by real-valued intervals). The coarse-grained representation allows us to extend a previously proposed algorithm for the Discretizable DGP (DDGP), the branch-and-prune (BP) algorithm. In the standard BP, atomic positions are fixed to unique positions at every node of the search tree: we rather represent atomic positions by a pair consisting of a feasible region, together with a most-likely position for the atom in this region. While the feasible region is a constant during the search, the associated position can be refined by considering the new distance constraints that appear at ...

## R029. 5th International Conference on Combinatorial Optimization and Applications, COCOA 2011

- Year: 2011
- Venue: Lecture Notes in Computer Science (including subseries Lecture Notes in Artificial Intelligence and Lecture Notes in Bioinformatics)
- DOI: n/a
- Cited by: 0
- Document type: Conference review
- Query hits: s1-core-ddgp, s3-symmetry-counts
- Tags hint: DDGP/DMDGP definitions, solution counting, distance geometry applications
- Authors:
- Author keywords: n/a
- Index keywords: n/a

Abstract: The proceedings contain 43 papers. The special focus in this conference is on Combinatorial Optimization and Applications. The topics include: Efficient algorithms for finding the k most vital edges for the minimum spanning tree problem; euclidean chains and their shortcuts; list dynamic coloring of sparse graphs; further improvement on maximum independent set in degree-4 graphs; approximation algorithms for minimum energy multicast routing with reception cost in wireless sensor networks; public communication based on russian cards protocol: A case study; minimum latency data aggregation in wireless sensor network with directional antenna; A near-optimal memoryless online algorithm for FIFO buffering two packet classes; on the maximum locally clustered subgraph and some related problems; algorithms for testing monomials in multivariate polynomials; quickest paths in anisotropic media; mechanisms for obnoxious facility game on a path; algorithmic aspects of heterogeneous biological networks comparison; minimum interval cover and its application to genome sequencing; exponential and ...

## R030. Clifford Algebra and the Discretizable Molecular Distance Geometry Problem

- Year: 2015
- Venue: Advances in Applied Clifford Algebras
- DOI: 10.1007/s00006-015-0532-2
- Cited by: 36
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications, algebraic/genericity methods
- Authors: Lavor C.; Alves R.; Figueiredo W.; Petraglia A.; Maculan N.
- Author keywords: distance geometry; Geometric algebra; protein conformation
- Index keywords: n/a

Abstract: Nuclear Magnetic Resonance experiments can provide distances between pairs of atoms of a protein that are close enough and the problem is how to determine the 3D protein structure based on this partial distance information, called Molecular Distance Geometry Problem. It is possible to define an atomic order 1,.., n and solve the problem iteratively using an exact method, called Branch-and-Prune (BP). The main step of BP algorithm is to solve a quadratic system to get the two possible positions for i, i > 3, in terms of the positions of i−3, i−2, i−1 and the distances di−1, i, di−2, i, di−3, i. Because of uncertainty in NMR data, some of the distances di−3, i may not be precise and the main problem to apply BP is related to the difficulty of obtaining an analytical expression of the position of atom i in terms of the positions of the three previous ones and the corresponding distances. We present such expression and although it is similar to one already existing in the literature, based on polyspherical coordinates, a new proof is given, based on Clifford algebra, and we also explain ...

## R031. Clifford algebra and discretizable distance geometry

- Year: 2018
- Venue: Mathematical Methods in the Applied Sciences
- DOI: 10.1002/mma.4422
- Cited by: 19
- Document type: Conference paper
- Query hits: s2-branch-prune, s5-lateration
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications, algebraic/genericity methods
- Authors: Alves R.; Lavor C.; Souza C.; Souza M.
- Author keywords: 3D protein structure; branch-and-prune algorithm; conformal Clifford algebra; distance geometry
- Index keywords: Algebra; Atoms; Bioinformatics; Geometry; Iterative methods; Proteins; Branch and prunes; Branch-and-prune algorithms; Clifford algebra; Combinatorial method; Distance geometry; Nuclear Magnetic Resonance (NMR); Protein structure calculation; Protein structures; Nuclear magnetic resonance

Abstract: Protein structure calculations using nuclear magnetic resonance (NMR) experiments are one of the most important applications of distance geometry. The chemistry of proteins and the NMR data allow us to define an atomic order, such that the distances related to the pairs of atoms {i−3,i},{i−2,i},{i−1,i} are available, and solve the problem iteratively using a combinatorial method, called branch-and-prune. The main step of BP algorithm is to intersect three spheres centered at the positions for atoms i−3,i−2,i, with radius given by the atomic distances di−3,i,di−2,i,di−1,i, respectively, to obtain the position for atom i. Because of uncertainty in NMR data, some of the distances di−3,i may not be precise or even not be available. Using conformal Clifford algebra, in addition to take care of NMR uncertainties, which implies that we have to calculate sphere intersections considering that their centers and radius may not be fixed anymore, we consider a more flexible atomic order, where distances di−3,i are replaced by dj,i, where j⩽i−3. Copyright © 2017 John Wiley & Sons, Ltd. Copyright ...

## R032. A least-squares approach for discretizable distance geometry problems with inexact distances

- Year: 2020
- Venue: Optimization Letters
- DOI: 10.1007/s11590-017-1225-7
- Cited by: 3
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications, algebraic/genericity methods
- Authors: Gonçalves D.S.
- Author keywords: Distance geometry; EDM; Gram matrix; Least-squares
- Index keywords: Additive noise; Geometry; Candidate positions; Distance geometry; Gram matrices; Least Square; Least squares solutions; Least-squares approach; Positive semidefinite; Quadratic equations; Least squares approximations

Abstract: A branch-and-prune (BP) algorithm is presented for the discretizable distance geometry problem in RK with inexact distances. The algorithm consists in a sequential buildup procedure where possible positions for each new point to be localized are computed by using distances to at least K previously placed reference points and solving a system of quadratic equations. Such a system is solved in a least-squares sense, by finding the best positive semidefinite rank K approximation for an induced Gram matrix. When only K references are available, a second candidate position is obtained by reflecting the least-squares solution through the hyperplane defined by the reference points. This leads to a search tree which is explored by BP, where infeasible branches are pruned on the basis of Schoenberg’s theorem. In order to study the influence of the noise level, numerical results on some instances with distances perturbed by a small additive noise are presented. © 2017, Springer-Verlag GmbH Germany, part of Springer Nature.

## R033. Integer Programming, Constraint Programming, and Hybrid Decomposition Approaches to Discretizable Distance Geometry Problems

- Year: 2022
- Venue: INFORMS Journal on Computing
- DOI: 10.1287/ijoc.2020.1039
- Cited by: 1
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications, algebraic/genericity methods
- Authors: MacNeil M.; Bodur M.
- Author keywords: Constraint programming; Decomposition algorithms; Discretization order; Distance geometry; Integer programming
- Index keywords: Binary trees; Combinatorial optimization; Constraint programming; Constraint theory; Geometry; Trees (mathematics); Undirected graphs; Branch-and-prune algorithms; Constraint programming; Decomposition algorithm; Discretization order; Discretizations; Distance geometry; Edge weights; Geometry problems; Integer Program- ming; Vertex ordering; Integer programming

Abstract: Given an integer dimension K and a simple, undirected graph G with positive edge weights, the Distance Geometry Problem (DGP) aims to find a realization function mapping each vertex to a coordinate in RK such that the distance between pairs of vertex coordinates is equal to the corresponding edge weights in G. The so-called discretization assumptions reduce the search space of the realization to a finite discrete one, which can be explored via the branch-and-prune (BP) algorithm. Given a discretization vertex order in G, the BP algorithm constructs a binary tree where the nodes at a layer provide all possible coordinates of the vertex corresponding to that layer. The focus of this paper is on finding optimal BP trees for a class of discretizable DGPs. More specifically, we aim to find a discretization vertex order in G that yields a BP tree with the least number of branches. We propose an integer programming formulation and three constraint programming formulations that all significantly outperform the state-of-the-art cutting-plane algorithm for this problem. Moreover, motivated by ...

## R034. Protein Loop Modeling via the Discretizable Distance Geometry Problem with Hydrogen-Based NMR Constraints

- Year: 2026
- Venue: ACS Omega
- DOI: 10.1021/acsomega.5c06422
- Cited by: 0
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Marques R.S.; Souza M.; Lavor C.
- Author keywords: n/a
- Index keywords: n/a

Abstract: Protein loop modeling remains a fundamental challenge in computational biology due to the inherent flexibility of loops and their critical role in biological functions. In this work, we employ a discrete distance geometry formulation, efficiently solved using the Branch-and-Prune algorithm, with a key innovation being the incorporation of hydrogen atoms into the model. Hydrogen atoms bonded to N and Cα in the protein backbone introduce additional geometric constraints, and their inclusion is particularly justified in the context of nuclear magnetic resonance (NMR) experiments, where short-range hydrogen–hydrogen distances can be detected and provide valuable structural information. By integrating these experimentally accessible constraints into the modeling process, we refine the representation of protein conformations. Computational experiments demonstrate that incorporating hydrogen atoms reduces the conformational space, leading to a more constrained and biologically realistic model. Comparisons with hydrogen-free formulations confirm that our approach improves agreement with ...

## R035. Improving the sampling process in the interval Branch-and-Prune algorithm for the discretizable molecular distance geometry problem

- Year: 2021
- Venue: Applied Mathematics and Computation
- DOI: 10.1016/j.amc.2020.125586
- Cited by: 5
- Document type: Article
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Lavor C.; Souza M.; Carvalho L.M.; Gonçalves D.S.; Mucherino A.
- Author keywords: Branch-and-Prune; Distance geometry; Protein structure
- Index keywords: Geometry; Nuclear magnetic resonance; Proteins; Branch-and-prune algorithms; Combinatorial method; Computational experiment; Mathematical details; Molecular distance geometry problem; Nuclear magnetic resonance(NMR); Protein molecules; Protein structures; Iterative methods

Abstract: Protein structure determination using Nuclear Magnetic Resonance (NMR) experiments is one of the most important applications of Distance Geometry, called the Molecular Distance Geometry Problem (MDGP). Using special atomic orders on the protein molecule, the MDGP can be solved iteratively using a combinatorial method, called Branch-and-Prune (BP). In order to deal with uncertainties of NMR data, there is an extension of the BP algorithm, called interval BP, where the idea is to sample values from the interval distances associated to such uncertainties. We propose a method to improve this sampling process, by reducing the interval of uncertain distances before taking the samples. All the mathematical details necessary to understand the proposal and its implementation are provided, along with some computational experiments that indicate the proposed strategy improves the interval BP algorithm. © 2020 Elsevier Inc.

## R036. Sphere intersection algorithms for Molecular Distance Geometry Problem

- Year: 2013
- Venue: Proceedings of the 2013 39th Latin American Computing Conference, CLEI 2013
- DOI: 10.1109/CLEI.2013.6670661
- Cited by: 1
- Document type: Conference paper
- Query hits: s2-branch-prune, s5-lateration
- Tags hint: Branch-and-Prune, distance geometry applications
- Authors: Santos C.; De Freitas R.; Salvatierra M.
- Author keywords: algorithms; branch-and-prune; euclidean distance; internal coordinates; linear systems; sphere intersection
- Index keywords: Algorithms; Linear systems; Three dimensional; Branch and prunes; Computational challenges; Computational experiment; Computational processing time; Euclidean distance; Internal coordinate systems; Internal coordinates; Molecular distance geometry problem; Spheres

Abstract: The problem of estimating the full three-dimensional structure of a molecule, determining the position in space of all the atoms that compose it, is called Molecular Distance Geometry Problem (MDGP). To do this from an incomplete set of distances is NP-hard computational problem, where to get a feasible solution in a reasonable execution time presenting interesting mathematical and computational challenges. In this work, continuous and discrete mathematical approaches to solve MDGP is revised, based on the analysis of two types of calculating of sphere intersection: solving nonlinear systems from interatomic Euclidean distance equations, or solving internal coordinate systems using matrix multiplication techniques. We adapted the Branch-and-Prune (BP) method considering four spheres intersection. Computational experiments using instances from PDB benchmark are performed, determining the 3D structure based on our theoretical assumptions in a competitive computational processing time. © 2013 IEEE.

## R037. An artificial backbone of hydrogens for finding the conformation of protein molecules

- Year: 2009
- Venue: Proceedings - 2009 IEEE International Conference on Bioinformatics and Biomedicine Workshops, BIBMW 2009
- DOI: 10.1109/BIBMW.2009.5332119
- Cited by: 9
- Document type: Conference paper
- Query hits: s1-core-ddgp, s2-branch-prune
- Tags hint: DDGP/DMDGP definitions, Branch-and-Prune, distance geometry applications
- Authors: Lavor C.; Mucherino A.; Liberti L.; Maculan N.
- Author keywords: n/a
- Index keywords: Bioinformatics; Combinatorial optimization; Molecules; Nuclear magnetic resonance spectroscopy; Branch and prunes; Combinatorial optimization problems; Exact algorithms; Molecular distance geometry problem; NMR experiments; Protein backbone; Protein molecules; SIMPLE method; Proteins

Abstract: NMR experiments can provide distances between pairs of hydrogens of a protein molecule. The problem of identifying the coordinates of such hydrogens by exploiting the information on the distances is a Molecular Distance Geometry Problem (MDGP). In a previous work, we defined an artificial backbone of hydrogens related to the protein backbones, where a particular ordering was given to the hydrogens. This ordering allows to formulate the MDGP as a combinatorial optimization problem, to which we refer as the Discretizable MDGP (DMDGP) and that we efficiently solve by an exact algorithm, the Branch and Prune (BP) algorithm. Once the coordinates of the hydrogens have been found, the problem of finding the remaining backbone atoms (N, C and C) is another MDGP. In this short paper, we propose a simple method for solving the MDGP related to the backbone atoms N, C and C of a protein, where the coordinates of the hydrogens previously found by the BP algorithm are exploited. ©2009 IEEE.

## R038. DDGP versus SMILE in newly diagnosed advanced natural killer/T-cell lymphoma: A randomized controlled, multicenter, open-label study in China

- Year: 2016
- Venue: Clinical Cancer Research
- DOI: 10.1158/1078-0432.CCR-16-0153
- Cited by: 137
- Document type: Article
- Query hits: s1-core-ddgp
- Tags hint: DDGP/DMDGP definitions
- Authors: Li X.; Cui Y.; Sun Z.; Zhang L.; Li L.; Wang X.; Wu J.; Fu X.; Ma W.; Zhang X.; Chang Y.; Nan F.; Li W.; Su L.; Wang J.; Xue H.; Zhang M.
- Author keywords: n/a
- Index keywords: Adolescent; Adult; Antineoplastic Combined Chemotherapy Protocols; Asparaginase; China; Cisplatin; Deoxycytidine; Dexamethasone; Disease-Free Survival; Etoposide; Female; Humans; Ifosfamide; Lymphoma, Extranodal NK-T-Cell; Male; Methotrexate; Middle Aged; Natural Killer T-Cells; Neoplasm Recurrence, Local; Polyethylene Glycols; Remission Induction; Young Adult; alanine aminotransferase; asparaginase; asparaginase macrogol; aspartate aminotransferase; cisplatin; creatinine; dexamethasone; etoposide; gemcitabine; ifosfamide; methotrexate; nitrogen; urea; antineoplastic agent; asparaginase; asparaginase macrogol; cisplatin; deoxycytidine; dexamethasone; etoposide; gemcitabine; ifosfamide; macrogol derivative; methotrexate; adult; advanced cancer; aged; allergy; alopecia; anemia; Article; cancer combination chemotherapy; China; clinical article; comparative effectiveness; controlled study; diarrhea; drug safety; drug tolerability; female; heart arrhythmia; heart failure; human; hyperbilirubinemia; hypertransaminasemia; hypofibrinogenemia; leukopenia; male; mucosa inflammation; multicenter study; multiple cycle treatment; nausea; neutropenia; NK T cell lymphoma; overall survival; partial thromboplastin time; priority journal; progression free survival; randomized controlled trial; thrombocytopenia; urea nitrogen blood level; vomiting; adolescent; analogs and derivatives; clinical trial; disease free survival; drug effects; middle aged; natural killer T cell; NK T cell lymphoma; procedures; remission; tumor recurrence; young adult

Abstract: Purpose: Optimal treatment strategies for advanced natural killer/T (NK/T)-cell lymphoma have not been fully defined. We compared the safety and efficacy of DDGP and SMILE regimens for advanced NK/T-cell lymphoma in a randomized controlled, multicenter, and open-label clinical trial. Experimental Design: Patients were newly diagnosed in stages III-IV and had performance scores in 0 to 2. Six cycles of DDGP (dexamethasone, cisplatin, gemcitabline, and pegaspargase) or SMILE (dexamethasone, methotrexate, ifosfamide, L-asparaginase, and etoposide) chemotherapy were randomly assigned to them. The primary end point was progression-free survival (PFS). Secondary end points included response rate and overall survival (OS). The trial is ongoing and is registered with ClinicalTrials.gov (No. NCT01501149). Results: Of 42 patients enrolled, 21 were treated with DDGP therapy, and 21 patients were treated with SMILE therapy. The 1-year PFS and 2-year OS rates were better in the DDGP group than that in the SMILE group (86% vs. 38% for 1-year PFS, P = 0.006; 74% vs. 45% for 2-year OS, P = 0.027). ...

## R039. A Branch-and-Prune algorithm for the Molecular Distance Geometry Problem

- Year: 2008
- Venue: International Transactions in Operational Research
- DOI: 10.1111/j.1475-3995.2007.00622.x
- Cited by: 115
- Document type: Article
- Query hits: s2-branch-prune
- Tags hint: Branch-and-Prune, distance geometry applications
- Authors: Liberti L.; Lavor C.; Maculan N.
- Author keywords: Branch-and-Prune algorithm; Molecular Distance Geometry Problem
- Index keywords: n/a

Abstract: The Molecular Distance Geometry Problem consists in finding the positions in of the atoms of a molecule, given some of the inter-atomic distances. We show that under an additional requirement on the given distances this can be transformed to a combinatorial problem. We propose a Branch-and-Prune algorithm for the solution of this problem and report on very promising computational results. © 2017 Wiley. All rights reserved.

## R040. Distance geometry: Theory, methods, and applications

- Year: 2013
- Venue: Distance Geometry: Theory, Methods, and Applications
- DOI: 10.1007/978-1-4614-5128-0
- Cited by: 76
- Document type: Book
- Query hits: s2-branch-prune
- Tags hint: Branch-and-Prune, partial reflection symmetries, distance geometry applications, algebraic/genericity methods
- Authors: Mucherino A.; Liberti L.; Lavor C.; MacUlan N.
- Author keywords: n/a
- Index keywords: Drawing (graphics); Genetic algorithms; Heuristic methods; Hyperbolic functions; Matrix algebra; Molecular biology; Molecular dynamics; Rigidity; Sensor networks; Simulated annealing; Surveys; Tabu search; Computational properties; Dimensional structures; Euclidean distance matrices; Molecular dynamics simulations; Rational drug designs; Semi-definite programming; Sensor network localization; Variable neighborhood search; Geometry

Abstract: This volume is a collection of research surveys on the Distance Geometry Problem (DGP) and its applications. It will be divided into three parts: Theory, Methods and Applications. Each part will contain at least one survey and several research papers. The first part, Theory, will deal with theoretical aspects of the DGP, including a new class of problems and the study of its complexities as well as the relation between DGP and other related topics, such as: distance matrix theory, Euclidean distance matrix completion problem, multispherical structure of distance matrices, distance geometry and geometric algebra, algebraic distance geometry theory, visualization of K-dimensional structures in the plane, graph rigidity, and theory of discretizable DGP: symmetry and complexity. The second part, Methods, will discuss mathematical and computational properties of methods developed to the problems considered in the first chapter including continuous methods (based on Gaussian and hyperbolic smoothing, difference of onvex functions, semidefinite programming, branch-and-bound), discrete ...
