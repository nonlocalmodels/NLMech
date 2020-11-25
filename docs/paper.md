---
title: 'NLMech: Implementation of finite difference/meshfree discretization of nonlocal fracture models'
tags:
  - Peridynamics
  - Finite difference
  - Finite element
  - HPX
  - Asynchronous many-task systems
authors:
  - name: Prashant K. Jha
    orcid: 0000-0003-2158-364X
    affiliation: 1
  - name: Patrick Diehl
    orcid: 0000-0003-0872-7098
    affiliation: 2
affiliations:
 - name: Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin, Austin, TX, United States of America
   index: 1
 - name: Center for Computation \& Technology, Louisiana State University, Baton Rouge, LA, United States of America
   index: 2
date: 25 November 2020
bibliography: paper.bib

---

![NLMech's logo which shows the obtained damage of a peridynamic simulation.\label{fig:logo}](../assets/logo/logo_joss.png)

# Summary

The open source code *NLMech* is an implementation of finite difference approximation of nonlocal models, \emph{e.g.}\ peridynamic. Peridynamic (PD) [@silling2007peridynamic,@silling2005meshfree] is a nonlocal formulation of classical continuum mechanics that is particularly robust with mechanical deformations involving crack (discontinuous displacement) and damage. The model seamlessly handles the two regime of deformation: elastic/smooth deformation and fracture. The constitutive laws describing the material behavior are simple to conseptualize and implement. Particularly, in numerical implementation, no special care for the modeling of cracks is required. Successful comparison of PD against a variety of experiments have been done [@diehl2019review]. 

Unlike classical continuum mechanics, where the internal force in the material is given by the stress at the material point, in PD, the internal force at given material point is due to the sum of the pairwise forces with the neighboring points. I.e. the force is expressed as the integral of the pairwise force density between the given point and another point in the neighborhood. Neighborhood of point $x$ is defined as set of nearby points $y$ sharing the bond with $x$. This is typically taken as ball centered at $x$ with radius $\delta$. $\delta$ is the nonlocal length scale and is referred to as \textit{horizon}. PD is often divided in two classes: bond-based and state-based models. In bond-based models, the two material points interact via a pairwise force law and the forces between the material points do not depend on the deformation state of surrounding points. In contrast, state-based models also describe the interaction between two material points but now the volumetric deformation in the neighborhood of two points play a role. The governing equation of motion for the bond-based PD [@silling2005meshfree] reads as

$$ \varrho(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) = \int\limits_{B_\delta(\mathbf{X})}\mathbf{f}(\mathbf{u}(t,\mathbf{X}')-\mathbf{u}(t,\mathbf{X}),\mathbf{X}'-\mathbf{X}) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D$$

and the governing equation for the state-based PD [@silling2007peridynamic] reads as 

$$  \varrho (\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) =  \int\limits_{B_\delta(\mathbf{X})} (T[\mathbf{X},t]\langle \mathbf{X}' - \mathbf{X} \rangle - T[\mathbf{X}',t]\langle \mathbf{X} - \mathbf{X}' \rangle) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D \text{.} $$
Here $\varrho$ denotes the material's density, $\mathbf{u}$ displacement field in the material, $\ddot{\mathbf{u}}$ acceleration, and $\mathbf{b}$ external force density. The constitutive law, relating bond strain with bond force, is prescribed using either the pairwise force function $\mathbf{f}$ or the PD state $T$ [@silling2007peridynamic]. In the library, following material models are implemented:

* Elastic state-based PD model [@silling2007peridynamic],
* Prototype micro-elastic brittle bond-based PD model [@silling2005meshfree],
* Nonlinear bond-based PD model [@lipton2014dynamic,@lipton2016cohesive], and
* Nonlocal double-well state-based peridynamic model [@Lipton2018,@jha2019numerical].

For example input files for these models, we refer to the collection of [examples](https://nonlocalmodels.github.io/examples/) in the documentation.

Currently, the library supports finite difference (or more generally meshfree) discretization. Using the triangulation of arbitrary domain using library such as Gmsh [@geuzaine2009gmsh], the library can create a meshfree discretization. Library is quipped with necessary modules, such as FE elements and quadrature integration rules, for finite element discretization of PD. We briefly discuss the finite difference/meshfree discretization of PD next. Figure \autoref{fig:discrete} shows the domain D discretized with the nodes $X = \{ X_i \in \mathbb{R}^3 \vert i=1,\ldots,n\}$. Each node $X_i$ represents a small area/volume denoted by $V_i$. In PD, as said earlier, each point $X_i$ interacts with neighboring points in discrete ball $B_\delta(X_i) = \{X_j: |X_i - X_j| < \delta \}$. 

![ Adpated from [@Diehl2020].\label{fig:discrete}](discrete.pdf)

The discrete equation of motion is written as, for the bond-based PD,

$$ \varrho(X_i)\ddot{\mathbf{u}}(t,X_i) = \sum\limits_{j \in B_\delta(X_i)}\mathbf{f}(\mathbf{u}(t,X_j)-\mathbf{u}(t,X_i),X_j-X_i) V_j + \mathbf{b}(t,X_i) \text{ in } D,$$

and, the state-based PD,

$$  \varrho (X_i)\ddot{\mathbf{u}}(t,X_i) =  \sum\limits_{j \in B_\delta(X_i)} (T[X_i,t]\langle X_j - X_i \rangle - T[X_j,t]\langle X_i - X_j \rangle) V_j + \mathbf{b}(t,X_i) \text{ in } D \text{.} $$

Here $\mathbf{u}(t,X_i)$ denotes the displacement of node $X_i$ at time $0 \leq t\leq T$. For the time discretization, we can consider: \textit{1)} implicit time integration and \textit{2)} explicit time integration using either central difference or velocity verlet scheme.

NLMech relies on the following open source softwares: HPX [@Kaiser2020], Blaze [@iglberger2012high], Blaze_Iterative, Gmsh [@geuzaine2009gmsh], VTK [@schroeder2004visualization], and yaml-cpp. For details 
about the specific version, we refer to NLMech's [documentation](https://github.com/nonlocalmodels/NLMech#building).

## Applications 

NLMech was used for the following applications/publications:

* Numerical convergence of finite difference approximations for state based perdidynamic fracture
models [@jha2019numerical] 
* Complex fracture nucleation and evolution with nonlocal elastodynamics [@lipton2019complex]
* Free damage propagation with memory [@lipton2018free] 
* Kinext relations and local energy balance for linear elastic fracture mechanics from a
nonlocal peridynamic model [@jha2020kinetic]

For an updated list of applications/publications, we refer to corresponding [NLMech documentation](https://nonlocalmodels.github.io/publications/).

# Statement of need

Nonlocal models, like peridynamic, are computational expensive, like molecular dynamics or smoothed-particle hydrodynamics. Several 
publications of GPU-based implementations [@mossaiby2017opencl,@diehl2012implementierung,@diehl2015efficient] and one commercial implementation in LS-DYNA [@ren20173d] can be found in literature. However, 
from an open source perspective only two other peridynamic implementations: [Peridigm](https://github.com/peridigm/peridigm) [@littlewood2015roadmap] and [PDLammps](https://lammps.sandia.gov/doc/pair_peri.html) [@parks2008implementing], are available. Both of these codes rely on the Message Passing Interface (MPI). On modern super computers' many core architectures where the threads per computational node increase, it is more and more important to focus on the fine-grain parallelism with increasing cores per computational nodes. NLMech utilizes the C++ standard library for parallelism and concurrency (HPX) [@Kaiser2020] to address this challenge. For more details about use utilization of asynchronous many-task systems, we refer to [@diehl2018implementation]. Second, the code implements the experimental nonlinear bond-based and state-based model. Nonlinear bond-based PD models [@lipton2016cohesive] enables one to show mathematical properties in rigorous way:

* limit of PD to classical fracture mechanics theory (LEFM) [@lipton2016cohesive] 
* well-posedness of PD solutions in spaces appropriate for numerical analysis and apriori error estimates of numerical discretization [@jha2018numerical,@jha2019numerical].

# Future directions

Towards massively parallel implementation of PD, the project titled "Domain decomposition, load balancing, and massively parallel solvers for the class of nonlocal models" was introduced in Google Summer of Code 2020. The goal of the project was to develop shared parallel library based on HPX for nonlocal diffusion equation and implement load balancing algorithm using novel ideas. With GSoC 2020 work as a base, we aim to extend NLMech from multithreading parallelism to massively parallel implementation. 

On the other hand, we are interested in implementation of new material models, new time discretization schemes, and local-nonlocal coupled formulation. 

# Acknowledgments

NLMech has been funded by:

*  Army Research Office Grant # W911NF-16-1-0456. PI was Dr. Robert Lipton at Louisiana State University. This grant supported Prashant K. Jha on a postdoctoral position from October 2016 - July 2019.
*  Canada Research Chairs Program under the Canada Research Chair in Multiscale Modelling of Advanced Aerospace Materials held by M. Lévesque; Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grants Program under Discovery Grant RGPIN-2016-06412.

For a updated list of previous and current funding, we refer to the corresponding [NLMech website](https://github.com/nonlocalmodels/NLMech#acknowledgements).

# References
