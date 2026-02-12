# Fully Discrete System
Now that we have described the spatial and temporal discretizations, we can put them together into a fully discrete system.

Euler Accelerate Half Step
$$
\begin{gather}
    \int_{K_i} e^n \phi \, \mathrm{d}V = \int_{K_i} \rho\mathcal{E}^n \phi \, \mathrm{d}V
        - \int_{K_i} \frac{1}{2\rho^n} \bm{p}^n \cdot \bm{p}^n \phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \bm{p}^{(1)} \phi \, \mathrm{d}V = \int_{K_i} \bm{p}^n\phi \, \mathrm{d}V
        + \frac{\Delta t}{2} \frac{q}{m} \int_{K_i} \rho^n \bm{E}^n\phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \mathcal{E}^{(1)} \phi \, \mathrm{d}V =
        \int_{K_i} \frac{1}{2 \rho^n} \bm{p}^{(1)} \cdot \bm{p}^{(1)} \phi \, \mathrm{d}V
        + \int_{K_i} e^n \phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
\end{gather}
$$

Euler Move Whole Step
$$
\begin{gather}
    \int_{K_i} \rho^{(2)} \phi \, \mathrm{d}V = \int_{K_i} \rho^{(1)} \phi \, \mathrm{d}V 
        - \Delta t \int_{\partial K_i} \mathcal{F}^{(1)}_{\rho} \phi \, \mathrm{d}A
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \bm{p}^{(2)} \phi \, \mathrm{d}V = \int_{K_i} \bm{p}^{(1)} \phi \, \mathrm{d}V
        - \Delta t \int_{\partial K_i} \mathcal{F}^{(1)}_{\bm{p}} \phi \, \mathrm{d}A
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \mathcal{E}^{(2)} \phi \, \mathrm{d}V = \int_{K_i} \mathcal{E}^{(1)} \phi \, \mathrm{d}V
        - \Delta t \int_{\partial K_i} \mathcal{F}^{(1)}_{\mathcal{E}} \phi \, \mathrm{d}A
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
\end{gather}
$$

Field Solve
$$
\begin{gather}
    \varepsilon \int_{\Omega} \nabla \Phi^{n+1} \cdot \nabla \phi \, \mathrm{d}V
        - \varepsilon \int_{\partial \Omega} \nabla \Phi^{n+1} \cdot \bm{n} \phi \, \mathrm{d}{A}
        = \frac{q}{m} \int_{\Omega} \rho^{(2)} \phi \, \mathrm{d} V
        \quad \forall \, \phi \in H^1_{\nabla}(\Omega)
\end{gather}
$$

Euler Accelerate Half Step
$$
\begin{gather}
    \int_{K_i} e^{(2)} \phi \, \mathrm{d}V = \int_{K_i} \rho\mathcal{E}^{(2)} \phi \, \mathrm{d}V
        - \int_{K_i} \frac{1}{2\rho^{(2)}} \bm{p}^{(2)} \cdot \bm{p}^{(2)} \phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \bm{p}^{n+1} \phi \, \mathrm{d}V = \int_{K_i} \bm{p}^{(2)}\phi \, \mathrm{d}V
        + \frac{\Delta t}{2} \frac{q}{m} \int_{K_i} \rho^{(2)} \bm{E}^{n+1}\phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
    \int_{K_i} \mathcal{E}^{n+1} \phi \, \mathrm{d}V =
        \int_{K_i} \frac{1}{2 \rho^{n+1}} \bm{p}^{n+1} \cdot \bm{p}^{n+1} \phi \, \mathrm{d}V
        + \int_{K_i} e^{(2)} \phi \, \mathrm{d}V
        \quad \forall \, \phi \in V^0(\Omega) \quad \forall \, K_i \in \Omega \\
\end{gather}
$$