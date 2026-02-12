# Euler Plasma Models

## Euler Maxwell Equations

The Euler-Maxwell equations for a single species are
$$
\begin{gather}
    \dot{\rho} + \nabla \cdot \bm{p} = 0 \\
    \dot{\bm{p}} + \nabla \cdot \left[\frac{1}{\rho} \bm{p} \otimes \bm{p} + P \textrm{I}\right] =
        \frac{q}{m} \left( \rho \bm{E} + \bm{p} \times \bm{B} \right) \\
    \dot{\mathcal{E}} + \nabla \cdot \left[ \left(\mathcal{E} + P\right) \frac{\bm{p}}{\rho} \right] =
        \frac{q}{m} \bm{p} \cdot \bm{E} \\
    \varepsilon \dfrac{\partial \bm{E}}{\partial t} - \mu^{-1} \nabla \times \bm{B} = - \bm{J} \\
    \dfrac{\partial \bm{B}}{\partial t} + \nabla \times \bm{E} = 0
\end{gather}
$$
where $\rho$ is the mass density, $\bm{p}$ the momentum density, $\mathcal{E}$ the total fluid energy density, $P$ the pressure,
$\bm{E}$ and $\bm{B}$ the electromagnetic fields, $q$ is the species charge, and $m$ is its mass.
The current density $\bm{J}$ is $\frac{q}{m} \bm{p}$.

## Euler Poisson Equations
An electrostatic version for the Euler plasma model is an Euler Poisson equation.

$$
\begin{gather}
    \dot{\rho} + \nabla \cdot \bm{p} = 0 \\
    \dot{\bm{p}} + \nabla \cdot \left[\frac{1}{\rho} \bm{p} \otimes \bm{p} + P \textrm{I}\right] =
        \frac{q}{m} \left( \rho \bm{E} + \bm{p} \times \bm{B} \right) \\
    \dot{\mathcal{E}} + \nabla \cdot \left[ \left(\mathcal{E} + P\right) \frac{\bm{p}}{\rho} \right] =
        \dfrac{q}{m} \bm{p} \cdot \bm{E} \\
    - \varepsilon \Delta \Phi = \frac{q}{m} \rho \\
    \bm{E} = - \nabla \Phi
\end{gather}
$$

where $\Phi$ is the electrostatic potential. The field $\bm{B}$ is constant.

## Primitive Variables
The equations above are written in conservative variables, but the primitive variables are needed sometimes.
They are related as follows
$$
\begin{gather}
    \rho = m n \\
    \bm{p} = \rho \bm{v} \\
    \mathcal{E} = \frac{1}{2 \rho} \bm{p} \cdot \bm{p} + e
\end{gather}
$$
where $n$ is the number density, $\bm{v}$ is the velocity, and $e$ is the internal energy density.
The pressure $P$ is related to the internal energy density through an equation of state.
The ideal gas equation of state gives the pressure as
$$
\begin{gather}
    P = (\gamma - 1) e,
\end{gather}
$$
where $\gamma$ is the adiabatic index.
A temperature can also found for an ideal gas,
$$
\begin{gather}
    T = \frac{P}{n k},
\end{gather}
$$
where $k$ is Boltzmann's constant.