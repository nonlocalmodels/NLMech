
<script type="text/x-mathjax-config">
    MathJax.Hub.Config({
      tex2jax: {
        skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
        inlineMath: [['$','$']]
      }
    });
  </script>
  <script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Brief introduction to the equations

The governing equations and the discretization is briefly introduced on this page. For more details we
refer to the following two references [1,2].

## Governing equations

Unlike classical continuum mechanics, where the internal force is written in terms of the stress, in PD, the internal force at a given material point is due to the sum of the pairwise forces of the neighboring points. i.e. the force in PD is expressed as the integral of the pairwise force density between the given point and another point in the neighborhood. Neighborhood of point `$x$` is typically taken as the ball of radius `$\delta$`, centered at $x$, where `$\delta$` is the nonlocal length scale and is referred to as horizon. PD is often divided in two classes: bond-based and state-based models. In bond-based models, the two material points interact via a pairwise force law and the forces between the material points do not depend on the deformation state of surrounding points. In contrast, in the state-based models the volumetric deformation in the neighborhood of two points plays a role in the pariwise force. The governing equation of motion for the bond-based PD [3] reads as

$$ \varrho(\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) = \int\limits_{B_\delta(\mathbf{X})}\mathbf{f}(\mathbf{u}(t,\mathbf{X}')-\mathbf{u}(t,\mathbf{X}),\mathbf{X}'-\mathbf{X}) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D$$

and the governing equation for the state-based PD [4] reads as 

$$  \varrho (\mathbf{X})\ddot{\mathbf{u}}(t,\mathbf{X}) =  \int\limits_{B_\delta(\mathbf{X})} (T[\mathbf{X},t]\langle \mathbf{X}' - \mathbf{X} \rangle - T[\mathbf{X}',t]\langle \mathbf{X} - \mathbf{X}' \rangle) d\mathbf{X}' + \mathbf{b}(t,\mathbf{X}) \text{ in } D \text{.} $$

Here `$\varrho$` denotes density of the material, `$\mathbf{u}$` displacement field in the material, `$\ddot{\mathbf{u}}$` acceleration, and `$\mathbf{b}$` external force density. The constitutive law, relating bond strain with bond force, is prescribed using either the pairwise force function $\mathbf{f}$ or the PD state $T$ [4].


## Discretization 

Currently, the library supports finite difference (or more generally meshfree) discretization. Using the triangulation of arbitrary domain, the library can create a meshfree discretization. The library is equipped with necessary modules, such as FE elements and quadrature integration rules, for finite element discretization of PD. Next, we briefly discuss the finite difference/meshfree discretization of PD. Figure 1 shows the domain $D$ discretized with the nodes `$X = \{ X_i \in \mathbb{R}^3 \vert i=1,\ldots,n\}$`. Each node $X_i$ represents a small area/volume denoted by `$V_i$`. In PD, as previously mentioned, each point `$X_i$` interacts with neighboring points in ball (discrete) `$B_\delta(X_i) = \{X_j: \vert X_i - X_j| < \delta \}$`. 



The discrete equation of motion is written as, for the bond-based PD,

$$ \varrho(X_i)\ddot{\mathbf{u}}(t,X_i) = \sum\limits_{j \in B_\delta(X_i)}\mathbf{f}(\mathbf{u}(t,X_j)-\mathbf{u}(t,X_i),X_j-X_i) V_j + \mathbf{b}(t,X_i) \text{ in } D,$$

and, the state-based PD,

$$  \varrho (X_i)\ddot{\mathbf{u}}(t,X_i) =  \sum\limits_{j \in B_\delta(X_i)} (T[X_i,t]\langle X_j - X_i \rangle - T[X_j,t]\langle X_i - X_j \rangle) V_j + \mathbf{b}(t,X_i) \text{ in } D \text{.} $$

Here `$\mathbf{u}(t,X_i)$` denotes the displacement of node `$X_i$` at time `$0 \leq t\leq T$`. For the time discretization, we can consider: \textit{1)} implicit time integration and \textit{2)} explicit time integration using either central difference or velocity verlet scheme.


## References

1. P. Diehl, P. K. Jha, H. Kaiser, R. Lipton, and M. Lévesque. An asynchronous and task-based implementation of peridynamics utilizing hpx—the C++ standard library for parallelism and concurrency. SN Applied Sciences, 2(12):2144, 2020, [10.1007/s42452-020-03784-x]({https://doi.org/10.1007/s42452-020-03784-x), [Preprint](https://arxiv.org/abs/1806.06917). 
2. Jha, PK, Lipton R. "Numerical convergence of finite difference approximations for state based peridynamic fracture models."  Computer Methods in Applied Mechanics and Engineering, 1 July 2019, 351(1), 184 - 225. [Link](https://doi.org/10.1016/j.cma.2019.03.024)
3. Silling, S. A., & Askari, E. (2005). A meshfree method based on the peridynamic model of164solid mechanics.Computers & Structures,83(17-18), 1526–1535. [https://doi.org/10.1651016/j.compstruc.2004.11.026](https://doi.org/10.1651016/j.compstruc.2004.11.026)
4. Silling, S. A., Epton, M., Weckner, O., Xu, J., & Askari. (2007). Peridynamic states and167constitutive modeling.Journal of Elasticity,88(2), 151–184.[https://doi.org/10.1007/168s10659-007-9125-1](https://doi.org/10.1007/168s10659-007-9125-1)