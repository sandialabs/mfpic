# Euler Poisson Time Discretization

Now that we have described the spatial discretization, the next step is to discretize this system in time.
We will describe the time discretization in strong form at first for clarity.

**QUESTION: Should the following be written in terms of velocity $\bm{v}$ instead of momentum $\bm{p}$ this would make the comparison with PIC more clear and discussing momentum hides internal energy discretization problems**

One thing to note is that the terms
$$
\begin{gather}
    \dot{\bm{p}} = \frac{q}{m} \rho \bm{E} \\
    -\varepsilon \Delta \Phi = \rho
\end{gather}
$$
describe an oscillator, known as the plasma oscillation.
This can be seen more directly by plugging in definition of $\bm{E}$ and taking a time derivative of the second equation.
$$
\begin{gather}
    \dot{\bm{p}} = \frac{q}{m} \rho \bm{E} \\
    \varepsilon \nabla \cdot \dot{\bm{E}} = \frac{q}{m} \dot{\rho} = -\frac{q}{m} \nabla \cdot \bm{p} \\
    \varepsilon \dot{\bm{E}} = - \frac{q}{m} \bm{p}
\end{gather}
$$
If we wish to control the energy exchange in this oscillation without excess growth or decay we must use a symplectic integrator.

The simplest symplectic integrator is the symplectic Euler time integrator.
$$
\begin{gather}
    \varepsilon \bm{E}^{n+1} = \varepsilon \bm{E}^{n} - \Delta t \frac{q}{m} \bm{p}^{n} \\
    \bm{p}^{n+1} = \bm{p}^n + \Delta t \frac{q}{m} \rho \bm{E}^{n+1}
\end{gather}
$$

In order to relate this back to the electrostatic potential formulation consider the time discretization of the continuity
equation.
If the continuity equation is discretized in time with forward Euler
$$
\begin{gather}
    \rho^{n+1} = \rho^n - \Delta t \nabla \cdot \bm{p}^n,
\end{gather}
$$
then we can take the divergence of the $\bm{E}$ evolution equation and substitute in the discretization of the continuity
equation.
$$
\begin{gather}
    \varepsilon \nabla \cdot \bm{E}^{n+1} = \varepsilon \nabla \cdot \bm{E}^{n}
        - \Delta t \frac{q}{m} \nabla \cdot \bm{p}^{n} \\
    \varepsilon \nabla \cdot \bm{E}^{n+1} = \varepsilon \nabla \cdot \bm{E}^{n}
        + \frac{q}{m} \left(\rho^{n+1} - \rho^n \right) \\
    -\varepsilon \nabla \cdot \nabla \Phi^{n+1} = -\varepsilon \nabla \cdot \nabla \Phi^n
        + \frac{q}{m} \left(\rho^{n+1} - \rho^n \right) \\
    -\varepsilon \Delta \Phi^{n+1} = \frac{q}{m} \rho^{n+1} - \varepsilon \Delta \Phi^n - \frac{q}{m} \rho^n \\
    -\varepsilon \Delta \Phi^{n+1} = \frac{q}{m} \rho^{n+1}
\end{gather}
$$
This used the assumption that we satisfied the potential formulation at time $n$, i.e.
$-\varepsilon \Delta \Phi^n = \frac{q}{m} \rho^n$.

Putting these together we get
$$
\begin{gather}
    \rho^{n+1} = \rho^n - \Delta t \nabla \cdot \bm{p}^n \\
    -\varepsilon \Delta \Phi^{n+1} = \frac{q}{m} \rho^{n+1} \\
    \bm{p}^{n+1} = \bm{p}^n + \Delta t \frac{q}{m} \rho \bm{E}^{n+1} \, .
\end{gather}
$$

Symplectic Euler can be symmetrized into a second order method using Strang splitting.
This second order is sometimes called leapfrog or Verlet.
In the electric field formulation this is written as
$$
\begin{gather}
    \bm{p}^{n+1/2} = \bm{p}^n + \frac{\Delta t}{2} \frac{q}{m} \rho \bm{E}^n \\
    \varepsilon \bm{E}^{n+1} = \varepsilon \bm{E}^{n} - \Delta t \frac{q}{m} \bm{p}^{n+1/2} \\
    \bm{p}^{n+1} = \bm{p}^{n+1/2} + \frac{\Delta t}{2} \frac{q}{m} \rho \bm{E}^{n+1} \, .
\end{gather}
$$
In the potential formulation with the continuity equation this second order time discretization is
$$
\begin{gather}
    \bm{p}^{n+1/2} = \bm{p}^n + \frac{\Delta t}{2} \frac{q}{m} \rho^n \bm{E}^n \\
    \rho^{n+1} = \rho^n - \Delta t \nabla \cdot \bm{p}^{n+1/2} \\
    -\varepsilon \Delta \Phi^{n+1} = \frac{q}{m} \rho^{n+1} \\
    \bm{p}^{n+1} = \bm{p}^{n+1/2} + \frac{\Delta t}{2} \frac{q}{m} \rho^{n+1} \bm{E}^{n+1} \, .
\end{gather}
$$
Note that this directly analogous to the traditional Verlet time integration scheme for electrostatic PIC simulations.
Verlet for electrostatic PIC is normally written as
$$
\begin{gather}
    \bm{v}^{n+1/2} = \bm{v}^n + \frac{\Delta t}{2} \frac{q}{m} \bm{E}^n \\
    \bm{x}^{n+1} = \bm{x}^n + \Delta t \bm{v}^{n+1/2} \\
    -\varepsilon \Delta \Phi^{n+1} = \frac{q}{m} \rho^{n+1} \\
    \bm{v}^{n+1} = \bm{v}^{n+1/2} + \frac{\Delta t}{2} \frac{q}{m} \bm{E}^{n+1} \, .
\end{gather}
$$
The particle position, $\bm{x}$, update is generally called the particle move and is a particle discretization of the continuity
equation.
The charge density $\frac{q}{m} \rho$ is computed from particle positions.
The particle velocity, $\bm{v}$, update is generally referred to as the particle accelerate.
To highlight the similarities with the PIC time discretization we will use the terms move and accelerate
to describe the different stages of this time discretization.

## Euler Accelerate
Let's take a closer look at the Euler momentum update, i.e. the Euler accelerate,
$$
\begin{gather}
    \bm{p}^{n+1/2} = \bm{p}^n + \frac{\Delta t}{2} \frac{q}{m} \rho^n \bm{E}^n \\
\end{gather}
$$
This corresponds to the fluid source terms in the original PDE.
$$
\begin{gather}
    \dot{\bm{p}} = \frac{q}{m} \rho \bm{E} + \frac{q}{m} \bm{p} \times \bm{B} \\
    \dot{\mathcal{E}} = \frac{q}{m} \bm{p} \cdot \bm{E}
\end{gather}
$$
For now we will drop the the $\bm{B}$ term, i.e. magnetic Lorentz force.
**Question/TODO: How do we do the equivalent of the Boris push for fluids. In order to step over the cyclotron frequency and maintain no work being done by the magnetic field.**
Note the source on the fluid energy is present because the Lorentz force is changing the kinetic energy of the fluid but it isn't
changing the internal energy of the fluid.
This can be seen as follows
$$
\begin{gather}
    \mathcal{E} = \frac{1}{2 \rho} \bm{p} \cdot \bm{p} + e \\
    \dot{\mathcal{E}} = \frac{1}{\rho} \bm{p} \cdot \dot{\bm{p}} + \dot{e} \\
    \frac{q}{m} \bm{p} \cdot \bm{E} = \frac{1}{\rho} \bm{p} \cdot \frac{q}{m} \rho \bm{E} + \dot{e} \\
    \frac{q}{m} \bm{p} \cdot \bm{E} = \frac{q}{m} \bm{p} \cdot \bm{E} + \dot{e} \\
    \dot{e} = 0 \, .
\end{gather}
$$
where $e$ is the internal energy density.
So this could be written as a modified system of ODEs,
$$
\begin{gather}
    \dot{\bm{p}} = \frac{q}{m} \rho \bm{E} \\
    \dot{e} = 0 \, .
\end{gather}
$$

### Problem with Forward Euler on Fluid Energy ODE
According to the Verlet/symplectic integrator described above the momentum should be updated with a forward Euler step.
One approach to updating the fluid energy might be apply forward Euler, i.e.
$$
\begin{gather}
    \bm{p}^{n+1} = \bm{p}^n + \Delta t \frac{q}{m} \rho \bm{E}^n \\
    \mathcal{E}^{n+1} = \mathcal{E}^n + \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n \, .
\end{gather}
$$

There is a problem with this approach however.
This approach does not keep the internal energy density constant.
$$
\begin{gather}
  e^{n+1} = \mathcal{E}^{n+1} - \frac{1}{2\rho} \bm{p}^{n+1} \cdot \bm{p}^{n+1} \\
  e^{n+1} = \mathcal{E}^n_{i} + \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n
    - \frac{1}{2\rho} \left(\bm{p}^n + \Delta t \frac{q}{m} \rho \bm{E}^n\right)
    \cdot \left(\bm{p}^n + \Delta t \frac{q}{m} \rho^n \bm{E}^n\right) \\
  e^{n+1} = e^n + \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n - \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n
    - \frac{\Delta t^2}{2\rho} \frac{q}{m} \rho \bm{E}^n \cdot \frac{q}{m} \rho \bm{E}^n \\
  e^{n+1} = e^n - \frac{\Delta t^2}{2} \frac{q^2}{m^2} \rho \bm{E}^n \cdot \bm{E}^n \, .
\end{gather}
$$
If
$$
\begin{gather}
    \frac{\Delta t^2}{2} \frac{q^2}{m^2} \rho \bm{E}^n \cdot \bm{E}^n > e^n,
\end{gather}
$$
then $e^{n+1} < 0$, which is an invalid state for the Euler equations.

### Option 1: Update Momentum and Postprocess Energy
A different approach would to compute the internal energy density at time $n$, update the momentum with forward Euler, and then
postprocess the fluid energy with the new kinetic energy in order to keep the internal energy density constant, i.e.
$$
\begin{gather}
    e^n = \mathcal{E}^n - \frac{1}{2\rho} \bm{p}^n \cdot \bm{p}^n \\
    \bm{p}^{n+1} = \bm{p}^n + \Delta t \frac{q}{m} \rho \bm{E}^n \\
    \mathcal{E}^{n+1} = \frac{1}{2 \rho} \bm{p}^{n+1} \cdot \bm{p}^{n+1} + e^n \, .
\end{gather}
$$
This would guarantee $e^{n+1} = e^n$.

This option relies on the fact that we aren't fluxing energy through the boundaries of a cell during this step.
If internal energy fluxed through the boundary and then this postprocessing was done, the internal energy fluxes would be
disregarded.

### Option 2: Modify Energy Update
A second option would be to modify the source term on the fluid energy equation.
$$
\begin{gather}
    \bm{p}^{n+1} = \bm{p}^n + \Delta t \frac{q}{m} \rho \bm{E}^n \\
    \mathcal{E}^{n+1} = \mathcal{E}^n + \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n
        + \frac{\Delta t^2}{2} \frac{q^2}{m^2} \rho \bm{E}^n \cdot \bm{E}^n \\
\end{gather}
$$

Checking the internal energy density,
$$
\begin{gather}
  e^{n+1} = \mathcal{E}^{n+1} - \frac{1}{2\rho} \bm{p}^{n+1} \cdot \bm{p}^{n+1} \\
  e^{n+1} = \mathcal{E}^n_{i} + \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n
    + \frac{\Delta t^2}{2} \frac{q^2}{m^2} \rho \bm{E}^n \cdot \bm{E}^n
    - \frac{1}{2\rho} \left(\bm{p}^n + \Delta t \frac{q}{m} \rho \bm{E}^n\right)
    \cdot \left(\bm{p}^n + \Delta t \frac{q}{m} \rho^n \bm{E}^n\right) \\
  e^{n+1} = e^n + \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n
    + \frac{\Delta t^2}{2} \frac{q^2}{m^2} \rho \bm{E}^n \cdot \bm{E}^n
    - \Delta t \frac{q}{m} \bm{p}^n \cdot \bm{E}^n
    - \frac{\Delta t^2}{2\rho} \frac{q}{m} \rho \bm{E}^n \cdot \frac{q}{m} \rho \bm{E}^n \\
  e^{n+1} = e^n
\end{gather}
$$
shows that the internal energy density remains constant.

This option is equivalent to discretizing
$$
\begin{gather}
    \dot{\mathcal{E}} = \frac{q}{m} \bm{p} \cdot \bm{E}
        + \frac{\Delta t}{2} \frac{q^2}{m^2} \rho \bm{E}^n \cdot \bm{E}^n \\
\end{gather}
$$
with forward Euler.
I don't like how this seems to encode the time discretization into the residual.
If you used a different time discretization the modification to the residual would need to change.

This option seems like it could possibly work in a monolithic time integration scheme, where we were fluxing fluid energy at the
same time as accelerating the fluid.

## Euler Move
Now consider the update of the continuity equation.
The Euler equations have two other divergence operators, which are traditionally updated simultaneously with the same time
discretization along with the continuity equation.
This means we should evolve them all together with the forward Euler method.
$$
\begin{gather}
    \rho^{n+1} = \rho^n - \Delta t \nabla \cdot \bm{p}^n \\
    \bm{p}^{n+1} = \bm{p}^n - \Delta t \nabla \cdot \left[\frac{1}{\rho^n} \bm{p}^n \otimes \bm{p}^n + P^n \textrm{I}\right] \\
    \mathcal{E}^{n+1} = \mathcal{E}^n
        - \Delta t \nabla \cdot \left[ \left(\mathcal{E}^n + P^n\right) \frac{\bm{p}^n}{\rho^n} \right] \\
\end{gather}
$$

It seems like we could substitute in a higher order time integration scheme here, e.g. SSPERK22.
We would want to check the symplectic nature of the time integration in this case, since the time discretization of the
continuity equation connects the potential formulation to the plasma oscillation.

## Complete Time Discretization
Putting the more complete description of the Euler accelerate and move into the Verlet time integrator gives the following time
integration scheme. Option 1 is used for the accelerate for now.

Euler Accelerate Half Step
$$
\begin{gather}
    e^n = \mathcal{E}^n - \frac{1}{2\rho^n} \bm{p}^n \cdot \bm{p}^n \\
    \bm{p}^{(1)} = \bm{p}^n + \frac{\Delta t}{2} \frac{q}{m} \rho^n \bm{E}^n \\
    \mathcal{E}^{(1)} = \frac{1}{2 \rho^n} \bm{p}^{(1)} \cdot \bm{p}^{(1)} + e^n
\end{gather}
$$

Euler Move Whole Step
$$
\begin{gather}
    \rho^{(2)} = \rho^{(1)} - \Delta t \nabla \cdot \bm{p}^{(1)} \\
    \bm{p}^{(2)} = \bm{p}^{(1)} - \Delta t \nabla \cdot \left[\frac{1}{\rho^{(1)}} \bm{p}^{(1)} \otimes \bm{p}^{(1)} + P^n \textrm{I}\right] \\
    \mathcal{E}^{(2)} = \mathcal{E}^{(1)}
        - \Delta t \nabla \cdot \left[ \left(\mathcal{E}^{(1)} + P^n\right) \frac{\bm{p}^{(1)}}{\rho^{(1)}} \right] \\
\end{gather}
$$

Field Solve
$$
\begin{gather}
    -\varepsilon \Delta \Phi^{n+1} = \frac{q}{m} \rho^{(2)} \\
    \bm{E}^{n+1} = - \nabla \Phi^{n+1}
\end{gather}
$$

Euler Accelerate Half Step
$$
\begin{gather}
    e^{(2)} = \mathcal{E}^{(2)} - \frac{1}{2\rho^{(2)}} \bm{p}^{(2)} \cdot \bm{p}^{(2)} \\
    \bm{p}^{n+1} = \bm{p}^{(2)} + \frac{\Delta t}{2} \frac{q}{m} \rho^{(2)} \bm{E}^{n+1} \\
    \mathcal{E}^{n+1} = \frac{1}{2 \rho^{(2)}} \bm{p}^{n+1} \cdot \bm{p}^{n+1} + e^{(2)}
\end{gather}
$$

Note that density, $\rho$, internal energy density, $e$, and pressure $P$ don't change durring the accelerate stage, so some of
the stages of these variables are equal.
That is
$$
\begin{gather}
    \rho^{(1)} = \rho^n \quad \rho^{(2)} = \rho^{n+1} \\
    e^{(1)} = e^n \quad e^{(2)} = e^{n+1} \\
    P^{(1)} = P^n \quad P^{(2)} = P^{n+1} \, .
\end{gather}
$$