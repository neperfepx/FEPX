.. _development_history:

Development History (1990-2020)
===============================

by Paul R. Dawson, June 2020.

Milestones in the Development of FEpX
-------------------------------------

Development of present code began in the late 1990's. The impetus was to incorporate elasticity in the existing viscoplastic constitutive framework. A number of attributes of the mechanical response can be addressed only if elasticity is part of the total constitutive description, and several of these were important to the research goals at that time. The inclusion of elasticity, however, fundamentally alters the computational approach because elasticity is based on changes in configuration of differential volumes, whereas plasticity requires knowledge only of current configurations. Further, the stiff mathematical character of the resulting system of equations necessitates greater care be exercised in integrating the equations over time and motivates the use of implicit schemes. While introducing added complexity in this regard, the inclusion of elasticity removes the constraint of incompressibility, which is particularly difficult to enforce within a robust framework and limits options available for the computational methodologies. Given these factors, the intent of the effort was to develop a code to support research investigations with the following defining specifications:

- incorporates anisotropic elasto-viscoplastic behaviors, especially in materials with low rate sensitivity of yielding  such as metals at low homologous temperature;

- embodies nonlinear kinematics, which are necessary for handling the large strains and large rotations inherent in plastic flow;

- utilizes state-based representation of properties, with attention to verifiable state descriptions at microscale as  motivated  by the inclusion of crystallographic texture (for yield surfaces) and evolution of texture (for strain induced anisotropy);

- exercises robust numerical methodologies, including implicit integration for stress,  element types capable of large strain deformations, and general boundary conditions;

- features a data parallel implementation with good scaling characteristics; and,

- is an expansible code framework to facilitate testing of alternative plasticity models and numerical methodologies.

Incorporation of a code with elasticity and nonlinear kinematics (which are tightly coupled) was  carried out principally by E. Marin [MARIN98a]_, [MARIN98b]_. The starting point for this effort was a code that utilized a viscoplastic model implemented in parallel with a hybrid finite element formulation developed by A. Beaudoin [BEAUDOIN95]_ (the hybrid formulation was an effective approach for dealing with the incompressibility constraint in the presence of plastic anisotropy). This code had been migrated to a version using Fortran and MPI and had 10-node tetrahedral elements available for robust simulations of large plastic strains [MIKA99]_. The approach taken by Marin followed a methodology developed for isotropic elasto-viscoplastic behavior developed about a decade earlier by Eggert [EGGERT88]_. Refinements were added later by N. Barton [BARTON99]_ for anisotropic elasticity.

Earlier Development Leading to the Present Version
--------------------------------------------------

The use of finite elements with crystal plasticity models can be organized into two broad categories that are defined by the relative sizes of grains and elements [DAWSON98]_.  One category is labeled, 'large scale', and is defined by the grains being much smaller than elements; the other is labeled, 'small scale', and is defined by the elements some part of a grain.  In large scale applications, an ensemble of grains underlies spatial points at the continuum scale and defines the properties of the continuum at that point. In small scale applications, the volume within a finite element is entirely of one grain and the material exhibits properties of a single crystal. The code development leading eventually to FEpX began with a large scale implementation in which crystal plasticity was embedded in an Eulerian, viscoplastic formulation was devised by K. Mathur [MATHUR89]_ to model steady-state metal forming applications. The potential capabilities of a code that incorporated crystal plasticity was pursued because it permitted computing the  evolving anisotropy associated with plastic yielding directly derived from crystallographic texture. Subsequent application of the approach for metal forming and geologic flow demonstrated that polycrystal models were viable for flow fields that could be idealized as steady and two-dimensional [MATHUR90]_ [KUMAR95]_ [DAWSON00]_.

The desire to model transient processes, such as sheet forming, motivated the major effort to develop a data parallel code. The code developed in this effort employed a proprietary version of Fortran [MATHUR95]_ that managed interprocessor communications, and enabled the simulation of fully three-dimensional forming processes [BEAUDOIN93]_ [BEAUDOIN94]_. An interest in applying the approach to small scale as well as large scale problems subsequently led to development of the hybrid finite element formulation for polycrystals [BEAUDOIN95]_. This milestone solidified the role of the finite element approach for investigating the role of grain interactions in polycrystal deformations and opened the door to investigating the strengths and limitations of various mean field assumptions (e.g. Taylor, Sachs, Relaxed Constraints, and Constrained Hybrid). Concurrently, the parallel computing landscape was rapidly evolving, and to take advantage of the introduction of new platforms, the code architecture was re-structured to employ Message Passing Interface (MPI) routines to conduct interprocessor communications. This Fortran/MPI version remained limited to purely viscoplastic behaviors, but was exploited to study texture evolution in polycrystals as well as development of intragrain misorientation distributions. This code was the starting point for development of FEpX and a re-focusing of the simulation priorities on small scale applications.

Application-Driven Expansion of Capabilities
--------------------------------------------

FEpX development in the decade following the launch of the first version of the present code was centered on support of investigations related to mechanical behaviors of polycrystals. Improvements were made in numerical procedures to improve robustness (namely, the nonlinear solver, quadrature rules,  and integration of state variables). Modifications were implemented to provide options in the loading histories that enabled better replication of experimental loading protocols. In particular, options to invoke more complex sequences of loading, unloading and reloading used in *in situ* loading, x-ray and neutron diffraction experiments were implemented. Capabilities for cyclic and multiaxial loading were added.

The extraction of data  related to the orientation of the crystallographic lattice from the simulation was of paramount importance. Routines were implemented to identify elements of the mesh whose lattice orientations lie near crystallographic fibers, a process referred to as 'light-up' in analogy to diffraction measurements. The output of FEpX was coordinated with a number of ancillary capabilities for manipulation of orientation-dependent variables (ODFPF), representation of anisotropic yield surfaces,  and execution of a virtual diffractometer.

Definition of the virtual polycrystals simulated with FEpX was initially limited to regular tessellations comprised of dodecahedral grains. Every grain was discretized with tetrahedral elements, typically numbering from 48 for a coarse representation to 1536 for more finely resolved grains. The use of other regular tessellations (cubic and truncated octahedral, in particular) were also explored [RITZ09]_. The coupling of FEpX with Neper greatly improved the representation of virtual polycrystals by allowing for irregular  Voronoi or Laguerre tesselations and facilitating re-meshing in simulations taken to large plastic strains [QUEY11]_.

Individuals contributing to these improvements include: D. Boyce, R. Carson, T. Han, M. Kasemer, T. Marin, A. Poshadel, R. Quey, and S.-L. Wong.

Source Sharing and Documentation
--------------------------------

By 2010 the use of FEpX over a decade in a variety of research projects motivated a push for standardization, version control, and sharing best done within a collaborative platform. A repository was established in 2012 together with documentation (users manual compiled by A. Mitch) for the needed input and possible outputs for FEpX. In concert with the establishing the code repository, numerous improvements were made in the organization of input and output data. A full description of the underlying theory and finite element implementation was posted on arXiv in 2015 [DAWSON15]_. Individuals contributing to this effort include A. Poshadel and M. Kasemer.

Extensions of FEpX
------------------

One of the specifications of FEpX was to provide an expansible code framework to facilitate testing of alternative plasticity models and numerical methodologies. Such efforts typically require substantial alterations to the code and are not intended to result in permanent changes to the baseline code. Examples of investigations of this nature include: a kinematic model with slip gradients [GERKEN08]_;  a continuous intragrain lattice representation [CARSON19]_, and a kinematic framework for twinning [KASEMER20]_.

| Dr. Paul R. Dawson
| Joseph C. Ford Professor of Engineering Emeritus
| Sibley School of Mechanical and Aerospace Engineering
| Cornell University

References
----------

.. [BARTON99] N. R. Barton, P. R. Dawson, and M. P. Miller. Yield strength asymmetry predictions from polycrystal plasticity. *Journal of Engineering Materials and Technology*, 121:230-239, 1999.

.. [BEAUDOIN93] A. J. Beaudoin, K. K. Mathur, P. R. Dawson, and G.C. Johnson. Three-dimensional deformation process simulation with explicit use of polycrystalline plasticity models. *International Journal of Plasticity*, 9:833-860, 1993.

.. [BEAUDOIN94] A. J. Beaudoin, P. R. Dawson, K. K. Mathur, U. F. Kocks, and D. A. Korzekwa. Application of polycrystal plasticity to sheet forming. *Computer Methods in Applied Mechanics and Engineering*, 117:49-70, 1994.

.. [BEAUDOIN95] A. J. Beaudoin, P. R. Dawson, K. K. Mathur, and U. F. Kocks. A hybrid finite element formulation for polycrystal plasticity with consideration of macrostructural and microstructural linking. *International Journal of Plasticity*, 11:501-521, 1995.

.. [CARSON19] R. A. Carson and P. R. Dawson. Formulation and Characterization of a Continuous Crystal Lattice Orientation Finite Element Method (LOFEM) and its Application to Dislocation Fields. *Journal of the Physics and Mechanics of Solids*, 126:1-26, 2019.

.. [DAWSON98] P. R. Dawson and E. B. Marin. Computational mechanics for metal deformation processes using polycrystal plasticity. In Erik van der Giessen and Theodore Y. Wu, editors, *Advances in Applied Mechanics*, 34:78-169. Academic Press, 1998.

.. [DAWSON00] P. R. Dawson and H.-R. Wenk. Texturing the upper mantle during convection. *Philosophical Magazine A*, 80(3):573-598, 2000.

.. [DAWSON15] P. R. Dawson and D. E. Boyce. FEpX – Finite Element Polycrystals: Theory, finite element formulation, numerical implementation and illustrative examples. *arXiv:1504.03296 [cond-mat.mtrl-sci]*, 2015.

.. [EGGERT88] G. M. Eggert and P. R. Dawson. A viscoplastic formulation with elasticity for transient metal forming. *Computer Methods in Applied Mechanics and Engineering*, 70:165-190, 1988.

.. [GERKEN08] J. M. Gerken and P. R. Dawson. A finite element formulation to solve a non-local constitutive model with stresses and strains due to slip gradients. *Computer Methods in Applied Mechanics and Engineering*, 197:1343-1361, 2008.

.. [KASEMER20] M. P. Kasemer and P. R. Dawson. A finite element methodology to incorporate kinematic activation of discrete deformation twins in a crystal plasticity framework. *Computer Methods in Applied Mechanics and Engineering*, 358:112653, 2020.

.. [KUMAR95] A. Kumar and P. R. Dawson. Polycrystal plasticity modeling of bulk forming with finite elements over orientation space. *Computational Mechanics*, 17:10-25, 1995.

.. [MARIN98a] E. B. Marin and P. R. Dawson. On modeling the elasto-viscoplastic response of metals using polycrystal plasticity. *Computer Methods in Applied Mechanics and Engineering*, 165:1-21, 1998.

.. [MARIN98b] E. B. Marin and P. R. Dawson. Elastoplastic finite element analysis of metal deformations using polycrystal constitutive models. *Computer Methods in Applied Mechanics and Engineering*, 165:23-41, 1998.

.. [MATHUR89] K. K. Mathur and P. R. Dawson. On modeling the development of crystallographic texture in bulk forming processes. *International Journal of Plasticity*, 5:67-94, 1989.

.. [MATHUR90] K. K. Mathur and P. R. Dawson. Texture development during wire drawing. *Journal of Engineering Materials and Technology*, 112(3):292-297, 1990.

.. [MATHUR95] K. K. Mathur. Parallel algorithms for large scale simulations in materials processing. In S. F. Shen and P. R. Dawson, editors, *Simulation of Materials Processing: Theory, Methods and Applications – NUMIFORM 95*, 109-114. A. A. Balkema, 1995.

.. [MIKA99] D. P. Mika and P. R. Dawson. Polycrystal Plasticity Modeling of Intracrystalline Boundary Textures. *Acta Materialia*, 47(4):1355-1369, 1999.

.. [QUEY11] R. Quey, P. R. Dawson, and F. Barbe. Large-scale 3-D random polycrystals for the finite element method: generation meshing and remeshing. *Computer Methods in Applied Mechanics and Engineering*, 200:1729-1745, 2011.

.. [RITZ09] H. Ritz and P. R. Dawson. Sensitivity to grain discretization of the simulated crystal stress distributions in fcc polycrystals. *Modeling and Simulation in Materials Science and Engineering*, 17:1-21, 2009.
