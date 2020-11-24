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
  - name: Prashant K. Jha
    orcid: 0000-0003-2158-364X
    affiliation: 2
affiliations:
 - name: Center for Computation \& Technology, Lousiana State University, Baton Rouge, LA, United States of America
   index: 1
 - name: Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin, Austin, TX, United States of America
   index: 2
date: 13 August 2017
bibliography: paper.bib

---

![NLMech's logo which shows the obtained damage of a peridynamic simulation.\label{fig:logo}](../assets/logo/logo_joss.png)

# Summary

The open source code *NLMech* is an implementation of finite element and finite difference approximation of nonlocal models, \emph{e.g.}\ peridynamic. Peridyanmic (PD) [@silling2007peridynamic,@silling2005meshfree] is a nonlocal formulation of classical continnum mechanics wit a focus on 
disconintuties as they arise at crack and fractures. Sucessfull comparision against a variaty of experiments have been done [@diehl2019review]. The governing equation of montion for bond-based PD [@silling2005meshfree] reads as


and the governing equaiton for state-based PD [@silling2007peridynamic] reads as 


The constitutive law of the material is described in the pair-wise force function f or the PD state T. Following materials models are implemented:

* Elastic state-based PD model [@silling2007peridynamic],
* Prototype micro-elastic brittle bond-based PD model [@silling2005meshfree],
* Nonlinear bond-based PD model [@lipton2014dynamic,@lipton2016cohesive], and
* Nonlocal Double Well state-based peridynamic model [@Lipton2018].

For example input files for these models, we refer to the collection of [examples](https://nonlocalmodels.github.io/examples/) in the documentation.

For the discretization in space two discretization schemes: 1) a finite difference approximation and 2) a finite element approximaiton is implemented. We briefly introduce the finite difference scheme and for the more complex finite element, we refer for the sake of keeping the paper short to [].


For the discretization in time following integraiton schemes are available: 

NLMech utilizes relies on following open source software: HPX [@Kaiser2020], Blaze [@iglberger2012high], Blaze_Iterative, Gmsh [@geuzaine2009gmsh], VTK [@schroeder2004visualization], and yaml-cpp. For details 
about the specific verison, we refer to NLMech's [documentation](https://github.com/nonlocalmodels/NLMech#building).

# Statement of need

Nonlocal models, like peridyanmics, are computaitonal expensive, like molecular dynamics or smoothed-particle hydrodynamics. Serveral 
publications of GPU-based implementations [@mossaiby2017opencl,@diehl2012implementierung,@diehl2015efficient] and one commerical implementation in LS-DYNA [@ren20173d] can be found in literature. However, 
from an open source perspective only two other peridynamic implementations: [Peridigm](https://github.com/peridigm/peridigm) [@littlewood2015roadmap] and [PDLammps](https://lammps.sandia.gov/doc/pair_peri.html) [@parks2008implementing], are available. Both of these codes rely on the Message Passing Interface (MPI). On modern super computers' many core archtitecutres where the threads per computaitonal node increase, it is more and more important to focus on the fine-grain parallism with increasing cores per computaitonal nodes. NLMech utilzes the C++ standard library for parallism and concurrency (HPX) [@Kaiser2020] to address this challenge. For more details about use utilzation of asynchronous many-task systems, we refer to [].
Second, the code implements the nonlinear bond-based and the nonlocal Double Well state-based model. The benefit of these models is that results from numerical analysis [@jha2019numerical,@jha2020kinetic] are available. Todo: Prashant elaborate here.

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
