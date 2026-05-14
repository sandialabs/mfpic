## Euler Poisson Equations
One electrostatic version for the Euler plasma model is given by,

$$
\begin{gather}
    \dot{\rho} + \nabla \cdot \bm{p} = 0 \\
    \dot{\bm{p}} + \nabla \cdot \left[\frac{1}{\rho} \bm{p} \otimes \bm{p} + P \textrm{I}\right] =
        \frac{q}{m} \rho \bm{E} \\
    \dot{\mathcal{E}} + \nabla \cdot \left[ \left(\mathcal{E} + P\right) \frac{\bm{p}}{\rho} \right] =
        \dfrac{q}{m} \bm{p} \cdot \bm{E} \\
    \varepsilon \dot{\bm{E}} = -\frac{q}{m} \bm{p} \\
\end{gather}
$$
Ampere's Law is often discared in favor of Gauss' Law which results in an equivalent system,
$$
\begin{gather}
    E = -\nabla \Phi \\
    -\varepsilon \Delta \Phi = \frac{q}{m} \rho
\end{gather}
$$

For these notes we will use the Ampere's Law version, as we are focusing on the plasma oscillation and don't require an
electrostatic potential.

## Operator Splitting
If we want to exactly conserve energy in the plasma oscillation we need to apply the Crank Nicolson time timestepping algorithm
to that mode.
In order to do that while still remaining explicit in the Euler solve, we must operator split the convection of the Euler fields
from the plasma oscillation.
The two operators are then the plasma oscillation,
$$
\begin{gather}
    \dot{\bm{p}} = \frac{q}{m} \rho \bm{E} \\
    \varepsilon \dot{\bm{E}} = -\frac{q}{m} \bm{p},
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
The Euler equations will be discretized in space with Discontinuous Galerkin, as described
[here](EulerPoissonSpatialDiscretization.md#euler-equations-in-dg) and will be integrated in time with Forward Euler as 
described [here](EulerPoissonTemporalDiscretization.md#euler-move).

## Plasma Oscillation Spatial and Temporal Discretization
The momentum density $\bm{p}$ is already discretized with Discontinuous Galerkin.
We would like $\bm{E}$ to be solved locally on a cell and easily couple with the momentum density, so we will use a Discontinuous
Galerkin representation for $\bm{E}$ as well.
Since both $\bm{p}$ and $\bm{E}$ are discretized the same in space it is possible to multiply by the DG test function and
integrate over a cell in order to be in weak form.
However we can multiply by the inverse of the cell mass matrix, so the ODEs can be solved directly on the degrees of the freedom
and the weak form is unnecessary.

In time we will use Crank-Nicolson in order to exactly conserved energy in the plasma oscillation,
$$
\begin{gather}
    \bm{p}^{n+1} = \bm{p}^n + \Delta t \frac{q}{m} \frac{\rho}{2} \left(\bm{E}^{n+1} + \bm{E}^n\right) \\
    \varepsilon \bm{E}^{n+1} = \varepsilon \bm{E}^n - \Delta t \frac{q}{m} \frac{1}{2} \left(\bm{p}^{n+1} + \bm{p}^n \right).
\end{gather}
$$
Rearranging so that $n+1$ are all on the LHS, gives
$$
\begin{gather}
    \bm{p}^{n+1} - \Delta t \frac{q}{m} \frac{\rho}{2} \bm{E}^{n+1} = \bm{p}^n + \Delta t \frac{q}{m} \frac{\rho}{2} \bm{E}^n \\
    \Delta t \frac{q}{m} \frac{1}{2} \bm{p}^{n+1} + \varepsilon \bm{E}^{n+1}
        = - \Delta t \frac{q}{m} \frac{1}{2} \bm{p}^n + \varepsilon \bm{E}^n.
\end{gather}
$$
If multiple species are present, so we have $\bm{p}_0$ and $\bm{p}_1$, then the system becomes
$$
\begin{gather}
    \bm{p}_0^{n+1} - \Delta t \frac{q_0}{m_0} \frac{\rho_0}{2} \bm{E}^{n+1}
        = \bm{p}_0^n + \Delta t \frac{q_0}{m_0} \frac{\rho_0}{2} \bm{E}^n \\
    \bm{p}_1^{n+1} - \Delta t \frac{q_1}{m_1} \frac{\rho_1}{2} \bm{E}^{n+1}
        = \bm{p}_1^n + \Delta t \frac{q_1}{m_1} \frac{\rho_1}{2} \bm{E}^n \\
    \Delta t \frac{q_0}{m_0} \frac{1}{2} \bm{p}_0^{n+1} + \Delta t \frac{q_1}{m_1} \frac{1}{2} \bm{p}_1^{n+1}
        + \varepsilon \bm{E}^{n+1}
        = - \Delta t \frac{q_0}{m_0} \frac{1}{2} \bm{p}_0^n - \Delta t \frac{q_1}{m_1} \frac{1}{2} \bm{p}_1^n
        + \varepsilon \bm{E}^n.
\end{gather}
$$

In matrix form
$$
\begin{gather}
    \begin{pmatrix}
        1                                    &        &                                              & - \Delta t \frac{q_0}{m_0} \frac{\rho_0}{2} \\
                                             & \ddots &                                              & \vdots \\
                                             &        & 1                                            & - \Delta t \frac{q_{s-1}}{m_{s-1}} \frac{\rho_{s-1}}{2} \\
        \Delta t \frac{q_0}{m_0} \frac{1}{2} & \cdots & \Delta t \frac{q_{s-1}}{m_{s-1}} \frac{1}{2} & \varepsilon \\
    \end{pmatrix}
    \begin{pmatrix}
        \bm{p}_0^{n+1} \\
        \vdots \\
        \bm{p}_{s-1}^{n+1} \\
        \bm{E}^{n+1}
    \end{pmatrix}
    =
    \begin{pmatrix}
        \bm{p}_0^n + \Delta t \frac{q_0}{m_0} \frac{\rho_0}{2} \bm{E}^n \\
        \vdots \\
        \bm{p}_1^n + \Delta t \frac{q_1}{m_1} \frac{\rho_1}{2} \bm{E}^n \\
        - \Delta t \frac{q_0}{m_0} \frac{1}{2} \bm{p}_0^n - \Delta t \frac{q_1}{m_1} \frac{1}{2} \bm{p}_1^n + \varepsilon \bm{E}^n
    \end{pmatrix}
\end{gather}
$$

