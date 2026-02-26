# Boltzmann Electrons
A common approximation in electrostatic PIC codes is the **adiabatic** or **Boltzmann electron** model.
Electrons are assumed to be in quasistatic equilibrium with the electrostatic potential $\varphi = \varphi(\vec{x})$
and therefore follow a Boltzmann distribution whereby
$$
f \propto \exp \left( \frac{e \varphi}{k T} \right),
$$
where $f = f(\vec{x}, \vec{v})$ is the electron single-particle distribution function,
$e$ is the elementary charge,
$k$ is the Boltzmann constant,
and $T$ is the temperature of the electrons.
Practically, a reference number density $n_0$ is used to normalize the spatial distribution $n = n(\vec{x})$ such that
$$
n = \int f(\vec{x}, \vec{v}) \, \mathrm{d} \vec{v} = n_0 \exp \left( \frac{e \varphi(\vec{x})}{k T} \right).
$$
Hence the charge density $\rho = \rho(\vec{x})$ is given by
$$
\rho = -e n_0 \left( \frac{e \varphi}{k T} \right),
$$
and Poisson's equation for the electrostatic potential becomes
$$
- \varepsilon_0 \Delta \varphi + e n_0 \left( \frac{e \varphi}{k T} \right) = \rho_i,
$$
where $\rho_i = \rho_i(\vec{x})$ is the charge density due to non-adiabatic species,
typically positively charged ions.

The weak form of this modified Poisson's equation is
$$
\epsilon_0 \sum_{j=0}^{N-1} \varphi_j \int_\Omega \nabla \psi_i \cdot \nabla \psi_i \, \mathrm{d} \vec{x}
- \epsilon_0 \int_{\partial \Omega} \psi_i \nabla \varphi \cdot \hat{n} \, \mathrm{d} S
+ e n_0 \int_\Omega \psi_i \exp \left( \frac{e}{k T} \sum_{j=0}^{N-1} \varphi_j \psi_j \right) \, \mathrm{d} \vec{x}
= \int_\Omega \psi_i \rho_i \, \mathrm{d} \vec{x}
$$
for basis functions $\psi_i$ of some finite element space $V_h \sub H^1(\Omega)$.
Unlike the usual Poisson's equation,
the stiffness operator is nonlinear.
The nonlinear form associated with the Boltzmann electron term is
$$
k(u, v) = e n_0 \int_\Omega u \exp \left( \frac{e v}{k T} \right) \, \mathrm{d} \vec{x},
$$
and its Jacobian is
$$
J(u,v) = \frac{e^2 n_0}{kT} \int_\Omega u \exp \left( \frac{e v}{k T} \right) \, \mathrm{d} \vec{x}.
$$

Poisson's equation with Boltzmann electrons is nonlinear and practically requires Newton's method to solve.
Given the finite element solution $\vec{\varphi}_n$ at Newton iteration $n$,
the linear system
$$
(L + J(\vec{\varphi}_n)) (\vec{\varphi}_{n+1} - \vec{\varphi}_n) = - \left( L \vec{\varphi}_n + \vec{k} (\vec{\varphi}_n)\right)
$$
is solved for $\vec{\varphi}_{n+1}$,
where $L$ is the discrete Laplacian matrix,
$$
k_i (\vec{u}) =  e n_0 \int_\Omega \psi_i \exp \left( \frac{e u_h}{k T} \right) \, \mathrm{d} \vec{x},
$$
$$
J_{ij}(\vec{u}) = \frac{e^2 n_0}{kT} \int_\Omega \psi_i \psi_j \exp \left( \frac{e u_h}{k T} \right) \, \mathrm{d} \vec{x},
$$,
and $u_h = \sum_{i=0}^{N-1} u_i \psi_i$.
