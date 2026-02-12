# Euler Poisson Spatial Discretization

In order to discretize the Euler Poisson equations we first describe two different finite element spaces.
The first is $H^1_{\nabla}(\Omega)$, where $\Omega$ is the domain of interest.
The domain $\Omega$ will be partitioned into a set of mesh elements, i.e. $\Omega = \left\{K_i\right\}$ where $K_i$ is mesh
element $i$.
The finite element space $H^1_{\nabla}(\Omega)$ is then defined as
$$
    H^1_{\nabla}(\Omega) = \left\{v \in C^0(\Omega) : v|_{K_i} \in P^1(K_i) \right\} \, .
$$
In words $H^1_{\nabla}(\Omega)$ is the set of functions that are continuous on $\Omega$ and linear on each element $K_i$.
Functions in the space have a well-defined gradient, hence why the space is labelled with $\nabla$.
The gradient may be discontinuous across element boundaries, but they will be integrable, that is
$$
    \nabla v \in L^2(\Omega) \quad \forall \quad v \in H^1_{\nabla}(\Omega) \, .
$$

The second finite element space is the discontinuous Galerkin, DG, finite element space, denoted $V^0(\Omega)$.
This space is defined as
$$
    V^0(\Omega) = \left\{v \in L^2(\Omega) : v|_{K_i} \in P^0(K_i) \right\} \, .
$$
In words $V^0$ is the set of integrable functions, that are piecewise constant on each element.
The DG finite element space is often modified into a vector finite element space, $\bm{V}^0(\Omega)$.
The space $\bm{V}^0$ is the space of vector functions where each component is in $V^0(\Omega)$, i.e.
$$
    \bm{V}^0(\Omega) = \left\{ \bm{v} : v_i \in V^0(\Omega) \right\} \, .
$$

To discretize the Euler Poisson equations in space we will let $\Phi \in H^1_{\nabla}(\Omega)$, $\rho \in V^0(\Omega)$,
$\rho \bm{v} \in \bm{V}^0(\Omega)$ and $\rho \mathcal{E} \in V^0(\Omega)$.
We now transform the partial differential equations above into weak form using these finite element spaces.

## Poisson Equation In HGrad

First consider the electrostatic potential equation sometimes called Gauss' Law.
The typical approach is to multiply by a test function in the finite element space, integrate over the domain, and then integrate
by parts where appropriate.
$$
\begin{gather}
    - \varepsilon \Delta \Phi = \frac{q}{m} \rho \\
    - \varepsilon \int_{\Omega} \Delta \Phi \phi \, \mathrm{d}V = \frac{q}{m} \int_{\Omega} \rho \phi \, \mathrm{d} V
        \quad \forall \, \phi \in H^1_{\nabla}(\Omega) \\
    \varepsilon \int_{\Omega} \nabla \Phi \cdot \nabla \phi \, \mathrm{d}V
        - \varepsilon \int_{\partial \Omega} \nabla \Phi \cdot \bm{n} \phi \, \mathrm{d}{A}
        = \frac{q}{m} \int_{\Omega} \rho \phi \, \mathrm{d} V
        \quad \forall \, \phi \in H^1_{\nabla}(\Omega) \, .
\end{gather}
$$
Note that because $\Phi \in H^1_{\nabla}(\Omega)$ the Laplacian is not well-defined on $\Phi$ - only the gradient is well-defined.
This is the reason why integration by parts is used to move the gradient onto the test function $\phi$.
This results in all differential operators being well-defined.
Also note that computing the right hand side is straightforward because $\rho$ piecewise constant.

Boundary conditions for this system must be considered in order to make this system solvable.
The "natural" boundary conditions for this problem are $\nabla \Phi \cdot \bm{n} = 0$ on $\partial \Omega$, this makes the
boundary term disappear, i.e.
$$
\begin{gather}
    \nabla \Phi \cdot \bm{n} |_{\partial \Omega} = 0 \Rightarrow
    \int_{\partial \Omega} \nabla \Phi \cdot \bm{n} \phi \, \mathrm{d}{A} = 0 \, .
\end{gather}
$$
If this is the only boundary condition used then infinitely many solutions can be generated from any given solution by adding a
constant to the potential.
In linear algebra terms the linear system becomes underdetermined.
In this case linear solver may fail.
To fix this issue numerically an additional constraint must be added, that effectively selects the constant.
The simplest constraint is known as pinning.
Pinning is performed by forcing $\Phi = 0$ for a specific node in the mesh.

The other boundary condition we consider for this system is a constant Dirichlet boundary condition, i.e.
$$
\begin{gather}
    \Phi |_{\Gamma} = c
\end{gather}
$$
where $\Gamma \subset \partial \Omega$.
Different constants can be used on different portions of the boundary and Dirichlet and natural boundary conditions can be mixed
and matched on different parts of the boundary as well.

## Euler Equations in DG

Now we consider the Euler equations with the DG formulation.
Note that normally the DG formulation is considered on each cell of the mesh individually.
First, let's consider the continuity equation
$$
\begin{gather}
    \dot{\rho} + \nabla \cdot \bm{p} = 0 \\
    \int_{K_i} \dot{\rho} \phi \, \mathrm{d}V + \int_{K_i} \nabla \cdot \bm{p} \phi \, \mathrm{d}{V} = 0
        \quad \forall \, \phi \in V^0(\Omega) \\
    \int_{K_i} \dot{\rho} \phi \, \mathrm{d}V - \int_{K_i} \bm{p} \cdot \nabla \phi \, \mathrm{d}{V}
        + \int_{\partial K_i} \widetilde{\bm{p} \cdot \bm{n}} \phi \, \mathrm{d}A = 0
        \quad \forall \, \phi \in V^0(\Omega) \\
    \int_{K_i} \dot{\rho} \phi \, \mathrm{d}V + \int_{\partial K_i} \mathcal{F}_{\rho} \phi \, \mathrm{d}A = 0
        \quad \forall \, \phi \in V^0(\Omega) \, .
\end{gather}
$$
Note that $\phi \in V^0(\Omega)$ and so is constant inside of an element and therefore the gradient is zero inside of each element
of the mesh.
Also note that $\bm{p}$ is discontinuous on element boundaries is therefore not well-defined on $\partial K_i$.
The tilde indicates that a numerical flux is required to evaluate the boundary term, the numerical flux will be denoted
$\mathcal{F}_*$ for different terms $*$ for conciseness.

We can do the same for the momentum and energy equations to get
$$
\begin{gather}
    \int_{K_i} \dot{\bm{p}} \phi \, \mathrm{d}V
        + \int_{\partial K_i} \mathcal{F}_{\bm{p}} \phi \, \mathrm{d}A =
        \dfrac{q}{m} \int_{K_i} \left(\rho \bm{E} + \bm{p} \times \bm{B} \right) \phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \\
    \int_{K_i} \dot{\mathcal{E}} \phi \, \mathrm{d}V
        + \int_{\partial K_i} \mathcal{F}_{\mathcal{E}} \phi \, \mathrm{d}A =
        \dfrac{q}{m} \int_{K_i} \bm{p} \cdot \bm{E} \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \, .
\end{gather}
$$

## Full Semi-Discrete System

All of these finite element equations can be put together to form a semi-discrete system.
$$
\begin{gather}
    \int_{K_i} \dot{\rho} \phi \, \mathrm{d}V + \int_{\partial K_i} \mathcal{F}_{\rho} \phi \, \mathrm{d}A = 0
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \dot{\bm{p}} \phi \, \mathrm{d}V
        + \int_{\partial K_i} \mathcal{F}_{\bm{p}} \phi \, \mathrm{d}A =
        \dfrac{q}{m} \int_{K_i} \left(\rho \bm{E} + \bm{p} \times \bm{B} \right) \phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \dot{\mathcal{E}} \phi \, \mathrm{d}V
        + \int_{\partial K_i} \mathcal{F}_{\mathcal{E}} \phi \, \mathrm{d}A =
        \dfrac{q}{m} \int_{K_i} \bm{p} \cdot \bm{E} \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \varepsilon \int_{\Omega} \nabla \Phi \cdot \nabla \phi \, \mathrm{d}V
        - \varepsilon \int_{\partial \Omega} \nabla \Phi \cdot \bm{n} \phi \, \mathrm{d}{A}
        = \frac{q}{m} \int_{\Omega} \rho \phi \, \mathrm{d} V
        \quad \forall \, \phi \in H^1_{\nabla}(\Omega)
\end{gather}
$$
