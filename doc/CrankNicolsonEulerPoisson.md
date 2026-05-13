## Euler Poisson Equations
An electrostatic version for the Euler plasma model is an Euler Poisson equation,

$$
\begin{gather}
    \dot{\rho} + \nabla \cdot \bm{p} = 0 \\
    \dot{\bm{p}} + \nabla \cdot \left[\frac{1}{\rho} \bm{p} \otimes \bm{p} + P \textrm{I}\right] =
        \frac{q}{m} \rho \bm{E} \\
    \dot{\mathcal{E}} + \nabla \cdot \left[ \left(\mathcal{E} + P\right) \frac{\bm{p}}{\rho} \right] =
        \dfrac{q}{m} \bm{p} \cdot \bm{E} \\
    - \varepsilon \Delta \Phi = \frac{q}{m} \rho \\
    \bm{E} = - \nabla \Phi
\end{gather}
$$

where $\Phi$ is the electrostatic potential.
In order to highlight the plasma oscillation, we can take the time derivative of the Poisson equation,
$$
\begin{gather}
    - \varepsilon \Delta \dot{\Phi} = \frac{q}{m} \dot{\rho} \\
    - \varepsilon \Delta \dot{\Phi} = -\frac{q}{m} \nabla \cdot \bm{p}.
\end{gather}
$$

## Operator Splitting
If we want to exactly conserve energy in the plasma oscillation we need to apply the Crank Nicolson time timestepping algorithm
to that mode.
In order to do that while still remaining explicit in the Euler solve, we must operator split the convection of the Euler fields
from the plasma oscillation.
The two operators are then the plasma oscillation,
$$
\begin{gather}
    \dot{\bm{p}} = -\frac{q}{m} \rho \nabla \Phi \\
    -\varepsilon \Delta \dot{\Phi} = -\frac{q}{m} \nabla \cdot \bm{p},
\end{gather}
$$
and the Euler move,
$$
\begin{gather}
    \dot{\rho} + \nabla \cdot \bm{p} = 0 \\
    \dot{\bm{p}} + \nabla \cdot \left[\frac{1}{\rho} \bm{p} \otimes \bm{p} + P \textrm{I}\right] = 0 \\
    \dot{\mathcal{E}} + \nabla \cdot \left[ \left(\mathcal{E} + P\right) \frac{\bm{p}}{\rho} \right] = 0.
\end{gather}
$$

## Crank Nicolson

$$
\begin{gather}
    -\frac{\varepsilon}{\Delta t} \Delta \left(\Phi^{n+1} - \Phi^n\right) =
        -\frac{q}{2 m} \nabla \cdot \left(\bm{p}^{n+1} + \bm{p}^n\right) \\
    \frac{1}{\Delta t} \left(\bm{p}^{n+1} - \bm{p}^n\right) = -\frac{q \rho}{2 m} \nabla \left(\Phi^{n+1} + \Phi^n\right)
\end{gather}
$$

$$
\begin{gather}
    -\varepsilon \Delta \Phi^{n+1} + \frac{q \Delta t}{2 m} \nabla \cdot \bm{p}^{n+1} = -\varepsilon \Delta \Phi^n
        - \frac{q \Delta t}{2 m} \nabla \cdot \bm{p}^n \\
    \frac{q \rho \Delta t}{2 m} \nabla \Phi^{n+1} + \bm{p}^{n+1} = - \frac{q \rho \Delta t}{2 m} \nabla \Phi^n + \bm{p}^n
\end{gather}
$$

Ampere Weak Form
$$
\begin{gather}
    - \varepsilon \int_{\Omega} \Delta \Phi^{n+1} \phi \, \mathrm{d}v
        + \frac{q \Delta t}{2 m} \int_{\Omega} \nabla \cdot \bm{p}^{n+1} \phi \, \mathrm{d}v =
        - \varepsilon \int_{\Omega} \Delta \Phi^n \phi \, \mathrm{d}v
        - \frac{q \Delta t}{2 m} \int_{\Omega} \nabla \cdot \bm{p}^n \phi \, \mathrm{d}v \\
    \varepsilon \int_{\Omega} \nabla \Phi^{n+1} \cdot \nabla \phi \, \mathrm{d}v
        - \varepsilon \int_{\partial \Omega} \nabla \Phi^{n+1} \cdot \bm{n} \phi \, \mathrm{d}a
        - \frac{q \Delta t}{2 m} \int_{\Omega} \bm{p}^{n+1} \cdot \nabla \phi \, \mathrm{d}v
        + \frac{q \Delta t}{2 m} \int_{\partial \Omega} \bm{p}^{n+1} \cdot \bm{n} \phi \, \mathrm{d}a =
        \varepsilon \int_{\Omega} \nabla \Phi^n \cdot \nabla \phi \, \mathrm{d}v
        - \varepsilon \int_{\partial \Omega} \nabla \Phi^n \cdot \bm{n} \phi \, \mathrm{d}a
        + \frac{q \Delta t}{2 m} \int_{\Omega} \bm{p}^n \cdot \nabla \phi \, \mathrm{d}v
        - \frac{q \Delta t}{2 m} \int_{\partial \Omega} \bm{p}^n \cdot \bm{n} \phi \, \mathrm{d}a \\
\end{gather}
$$
with periodic boundary conditions
$$
\begin{gather}
    \varepsilon \int_{\Omega} \nabla \Phi^{n+1} \cdot \nabla \phi \, \mathrm{d}v
        - \frac{q \Delta t}{2 m} \int_{\Omega} \bm{p}^{n+1} \cdot \nabla \phi \, \mathrm{d}v =
        \varepsilon \int_{\Omega} \nabla \Phi^n \cdot \nabla \phi \, \mathrm{d}v
        + \frac{q \Delta t}{2 m} \int_{\Omega} \bm{p}^n \cdot \nabla \phi \, \mathrm{d}v \quad \forall \quad \phi \in H_{\nabla}(\Omega)
\end{gather}
$$
Lorentz Weak Form
$$
\begin{gather}
    \frac{q \rho \Delta t}{2 m} \int_{\Omega} \nabla \Phi^{n+1} \varphi \, \mathrm{d}v
        + \int_{\Omega} \bm{p}^{n+1} \varphi \, \mathrm{d}v
        = \int_{\Omega} \bm{p}^n \varphi \, \mathrm{d}v
        - \frac{q \rho \Delta t}{2 m} \int_{\Omega} \nabla \Phi^n \varphi \, \mathrm{d}v
\end{gather}
$$