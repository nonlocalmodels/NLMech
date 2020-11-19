---
title: 'NLMech: Implementation of finite element and finite difference approximation of Nonlocal models'
tags:
  - Peridynamics
  - Finite difference
  - Finite element
  - HPX
  - Asynchronous many-task systems
authors:
  - name: Patrick Diehl
    orcid: 0000-0003-0872-7098
    affiliation: 1
affiliations:
 - name: Center for Computation \& Technology, Lousiana State University 
   index: 1
date: 13 August 2017
bibliography: paper.bib

---

![NLMech's logo which shows the obtained damage of a peridynamic simulation. For the input files and more details we refer to the collection of [examples](https://nonlocalmodels.github.io/examples/).\label{fig:logo}](../assets/logo/logo_sims_transparent.png)


# Summary


Following materials models are implemented:

* Elastic state-based PD model [@silling2007peridynamic],
* Prototype micro-elastic brittle bond-based PD model [@silling2005meshfree],
* Nonlinear bond-based PD model [@lipton2014dynamic,@lipton2016cohesive], and
* Nonlocal Double Well state-based peridynamic model [@Lipton2018].

For example input files for these models, we refer to the collection of [examples](https://nonlocalmodels.github.io/examples/) in the documentation.

# Statement of need

Nonlocal models, like peridyanmics, are computaitonal expensive, like molecular dynamics or smoothed-particle hydrodynamics. Serveral 
publications of GPU-based implementations [@mossaiby2017opencl,@diehl2012implementierung,@diehl2015efficient] and one commerical implementation in LS-DYNA [@ren20173d] can be found in literature. However, 
from an open source perspective only two other peridynamic implementations: [Peridigm](https://github.com/peridigm/peridigm) [@littlewood2015roadmap] and [PDLammps](https://lammps.sandia.gov/doc/pair_peri.html) [@parks2008implementing], are available. Both of these codes rely on the Message Passing Interface (MPI).

# Applications 

NLMech was used for following applications:

* Numerical convergence of finite difference approximations for state based perdidynamic fracture
models [@jha2019numerical] 
* Complex fracture nucleation and evoluation with nonlocal elastodynamics [@lipton2019complex]
* Free damage propogation with memory [@lipton2018free] 
* Kinext relations and local energy balance for linear elastic fracture mechancis from a
nonlocal peridynamic model [@jha2020kinetic]

For a updated list of applications, we refer to corresponding [NLMech documentation](https://nonlocalmodels.github.io/publications/).

# Acknowledgements

NLMech has been funded by:

*  Army Research Office Grant # W911NF-16-1-0456
*  Canada Research Chairs Program under the Canada Research Chair in Multiscale Modelling of Advanced Aerospace Materials held by M. Lévesque; Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grants Program under Discovery Grant RGPIN-2016-06412.

For a updated list of previous and current funding, we refer to the corresponding [NLMech website](https://github.com/nonlocalmodels/NLMech#acknowledgements).

# References
