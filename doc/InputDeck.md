# User Input

## Units

mfpic uses SI units. When appropriate, the expected units are listed below in brackets, `[-]`.

## Electrostatic Fields

Details for the electrostatic solve are specified under the ``Fields`` key.
The only option at the top level is the basis order given with key ``Basis Order``.

```yaml
Fields:
    Basis Order: int
```

### Electrostatic Boundary Conditions

Boundary conditions for the electrostatic solve are given under the ``Boundary Conditions`` key under ``Fields``.
Only constant dirichlet boundary conditions are available to be specified.
Each boundary condition should be specified as a item in a list, and give the side on which the boundary is to be applied and
the constant value of the boundary.
For inline meshes the sides can be specified as "left" $x = 0$, "right" $x = L_x$, "bottom" $y = 0$, "top" $y = L_y$,
"back" $z = 0$, or "front" $z = L_z$.
For meshes loaded from files the side is an integer specifying mesh attribute or sideset index.

If a boundary condition isn't specified on a side, then that side will have natural boundary conditions.
If all boundaries are natural then pinning will be applied to keep the problem from being underdetermined.

```yaml
Fields:
    Boundary Conditions:
        - Side: string or integer
          Value: double #[V]
```

**``Side``**: The side on which to apply the boundary condition.
    For inline meshes the side should be a string.
    The available sides are "left", "right", "top", "bottom", "front", and "back".
    For meshes loaded from files the side should be an integer.
    The integer is a mesh "attribute" in MFEM parlance or the index of the sideset in Cubit terminology.

**``Value``**: The constant value of the boundary condition.

#### Example

```yaml
Fields:
    Boundary Conditions:
        - Side: left
          Value: 1.
        - Side: right
          Value: 2.
```
specifies a problem with the potential fixed at 1V on the left boundary and 2V on the right boundary.

## Mesh

The mesh for the problem can be created in two ways.
A mesh can be created inline through user inputs or a mesh can be read in from a file.

### Inline Mesh

An inline mesh can be created from the input deck alone.
The meshes can be 1, 2, or 3 dimensional.
In 2D the domain is a rectangle given by $\Omega = [0, L_x] \times [0, L_y]$ and in 3D the domain is a box given by
$\Omega = [0, L_x] \times [0, L_y] \times [0, L_z]$.
The domain is then partitioned evenly in each dimension.
In 2D the rectangle is split into $N_x$ rectangles in the x direction and $N_y$ rectangles in the y direction.
If a `tri` mesh is specified then each rectangle is then split into two triangles.
In 3D the box is similarly partitioned with $N_z$ elements in the z direction as well.
If a `tet` mesh is specified then each box element is split into 6 tetrahedrons.

If a periodic mesh is desired, this can be set using the option ``Periodic Dimensions``.
This key is optional and can be omitted if periodicity is not required.


```yaml
Mesh:
    Type: string
    Lengths: list of doubles #[m]
    Number of Elements: list of integers
    Periodic Dimensions: list of strings
```

**``Type``**: The element type of the mesh, the available options are `line`, `quad`, `tri`, `hex`, and `tet`.

**``Lengths``**: The length of the mesh in each dimension.
    The meshes always start at 0.0 and extend to the length specified.

**``Number of Elements``**: The number of tensor mesh elements in each direction.
    If `tri` or `tet` is specified, then these elements will be split.
    For `tri` meshes each quadrilateral will be split into 2 triangles, and for `tet` meshes each hexahedron will be split into 6
    tetrahedrons.

**``Periodic Dimensions``**: A list of strings describing the dimensions that should be made periodic.
    The valid strings are `x`, `y`, and `z`.


#### Example
```yaml
Mesh:
    Type: quad
    Lengths: [0.2, 0.3]
    Number of Elements: [10, 15]
    Periodic Dimensions: [y]
```
This creates a square mesh, $\Omega = [0.0, 0.2] \times [0.0, 0.3]$, of quadrilaterals with $10$ elements in the x direction and
$15$ elements in the y direction.
This mesh is periodic in the y direction, that is the nodes at $y = 0.0$ and $y = 0.3$ are the same entities.

### File Mesh

A mesh can also be read in from a file.
The read in mesh can be made periodice through the option ``Periodic Dimensions``. Note that the sides in this direction must be
parallel and at the maximum length of the mesh in that dimension.

```yaml
Mesh:
    File Name: string
    Periodic Dimensions: list of strings
```

**``File Name``**: The name of the file to be loaded in as the mesh.

**``Periodic Dimensions``**: A list of strings describing the dimensions that should be made periodic.
    The valid strings are `x`, `y`, and `z`.

## Time Stepping

A constant timestep is required for simulations. The timestep size and number of timesteps can be determined by setting two of
the following three options, ``Final Time``, ``Number of Time Steps``, and ``Time Step Size``.
``Number of Time Steps`` and ``Time Step Size`` take precedence if all three are specified.
If ``Final Time`` and ``Time Step Size`` are specified then the number of timesteps and timestep size used in the simulation are
computed in a slightly complicated way.
The number of time steps is computed as
$$
    \text{Number of Time Steps} = \left\lceil\frac{\text{Final Time}}{\text{Time Step Size}}\right\rceil,
$$
and then the time step size actually used in the simulation is
$$
    \text{Time Step Size} = \frac{\text{Final Time}}{\text{Number of Time Steps}}.
$$

```yaml
Time Stepping:
    Final Time: double #[s]
    Number of Time Steps: integer
    Time Step Size: double #[s]
```

**``Final Time``**: The physical time that the simulation should end in seconds.

**``Number of Time Steps``**: The number of time steps to take in the simulation.

**``Time Step Size``**: The size of each time step in seconds.

Example
```yaml
Time Stepping:
    Number of Time Steps: 100
    Time Step Size: 0.12
```

## Species
Parameters for the species are added as maps under the `Species` key,
each key of which names the species internally.
```yaml
Species:
    Species Name:
        Mass: double #[kg]
        Charge: double #[C]
        Charge Over Mass: double (optional) #[C/kg]
        Specific Heat Ratio: double (default = 5/3)
```

**``Mass``**: Mass of each physical particle of the given species. Must be positive.

**``Charge``**: Charge of each physical particle of the given species.

**``Charge Over Mass``**: Override the charge-to-mass ratio used to accelerate members of this species.
If left unset, this is computed automatically from the given `Charge` and `Mass`.

**``Specific Heat Ratio``**: Specific heat ratio of the species.
Defaults to $5/3$, which is the specific heat ratio of a monatomic gas.

#### Example
```yaml
Species:
    electron:
        Mass: 9.11e-31
        Charge: -1.6e-19
    non-accelerated proton:
        Mass: 1.67e-27
        Charge: 1.6e-19
        Charge Over Mass: 0.0
```
This specifies two species, one of which is named `electron` and has roughly corresponding parameters,
and the other of which is named `non-accelerated proton` which contributes charge to the fields
but is not accelerated by them.

## Particles
Details for the (macro)particle physics are defined under the `Particles` key.

### Particle Initial Conditions
The initial particle populations are defined under the `Initial Conditions` key.
An arbitrary number of populations may be provided as a sequence,
each of which is a drifting Maxwellian distributed uniformly throughout the mesh.
Only constant-in-space initial conditions are presently supported.
```yaml
Particles:
    Initial Conditions:
        - Species: list of strings
          Number of Macroparticles per Species: integer
          Constant:
              Bulk Velocity: list of doubles, length 3 (default: [0.0, 0.0, 0.0]) #[m/s]
              Temperature: double #[K]
              Number Density: double #[-/m^3]
        ...
```

**``Species``**: List of species for which particles will be loaded.
These species must be defined in the top-level `Species` block.

**``Bulk Velocity``**: The bulk velocity of the particle distributions.
If left unset, this will be zero.

**``Temperature``**: Temperature of the particle distribution.
The temperature is computed independently for each species.
If set to zero, each particle will have its velocity set to `Bulk Velocity`.

**``Number Density``**: Number density of each species.

**``Number of Macroparticles per Species``**: Total number of particles loaded for each species.
The weight of each particle is given by
$$
    \text{Weight} = \frac{\text{Number Density} \times \text{Volume of Mesh}}{\text{Number of Macroparticles per Species}}.
$$

#### Example
```yaml
Particles:
    Initial Conditions:
        - Species: [electron]
          Number of Macroparticles per Species: 100
          Constant:
              Bulk Velocity: [1.0, 20.0, 300.0]
              Temperature: 0.0
              Number Density: 1.0e15
        - Species: [proton, electron]
          Number of Macroparticles per Species: 10000
          Constant:
              Temperature: 11600.0
              Number Density: 1.0e19
```

### Particle Sources
Particles can be dynamically loaded into the simulation with the optional `Sources` key.
Particles are sourced every time step.
As with the [Particle Initial Conditions](#particle-initial-conditions),
particles are loaded uniformly throughout the domain
according to a drifting Maxwellian velocity distribution.
The options are identical to those of the initial conditions;
the number density and number of macroparticles are loaded each time step.
Only constant-in-space sources are presently supported.

#### Example
```yaml
Particles:
    Sources:
        - Species: [proton, electron]
          Number of Macroparticles per Species: 10
          Constant:
              Temperature: 11600.0
              Number Density: 1.0e16
```

### Particle Boundary Conditions
The particle boundary conditions are defined under the `Boundary Conditions`
and `Default Boundary Condition` keys.
The value of `Boundary Conditions` must be a sequence,
each member of which specifies the type of boundary and the mesh side to which it applies.
The `Default Boundary Condition` is the type of boundary applied to every mesh side
for which a boundary condition has not explicitly been defined.
```yaml
Particles:
    Boundary Conditions:
        - Side: string or integer
          Type: boundary condition type string (options below)
        ...
    Default Boundary Condition: boundary condition type string (options below)
```

**``Side``**: The side on which to apply the boundary condition.
For inline meshes the side should be a string.
The available sides are "left", "right", "top", "bottom", "front", and "back".
For meshes loaded from files the side should be a positive integer.
The integer is a mesh "attribute" in MFEM parlance or the index of the sideset in Cubit terminology.

**``Type``**: The type of particle boundary.
The options are:
- **``Reflecting``**: Specularly reflects incident particles.
- **``Absorbing``**: Incident particles are instantly removed from the simulation.

Periodic boundaries are detected and applied automatically.

#### Example
```yaml
Particles:
    Boundary Conditions:
        - Side: left
          Type: Absorbing
        - Side: 3
          Type: Reflecting
    Default Boundary Condition: Reflecting
```

## Euler Fluids

Fluids that obey the Euler equations can be added under the ``Euler Fluids`` key.
The only option at the top level is ``Basis Order`` which specifies the discontinuous Galerkin finite element basis order that
will be used for all Euler fluids.

```yaml
Euler Fluids:
    Basis Order: int
```

### Euler Fluids Initial Conditions
The fluids present in the simulation and their initial conditions are specied under ``Initial Conditions``.
Each initial condition is an item in a list under ``Initial Conditions`` and can have more than one species.
Only constant-in-space or Sod shock tube initial conditions are supported.
If a constant in space is desired, then the ``Constant`` keyword and subparameters should be specified.
If a Sod shock tube problem is desired, then the ``Sod`` keyword and subparameters should be specified.
``Constant`` and ``Sod`` keys can not be used together.

```yaml
Euler Fluids:
    Initial Conditions:
        - Species: list of strings
          Constant:
            Bulk Velocity: list of doubles, length 3 (default: [0.0, 0.0, 0.0]) #[m/s]
            Temperature: double #[K]
            Number Density: double #[-/m^3]
          Sod:
            Discontinuity Location: double #[m]
            Left State:
                Bulk Velocity: list of doubles, length 3 (default: [0.0, 0.0, 0.0]) #[m/s]
                Temperature: double #[K]
                Number Density: double #[-/m^3]
            Right State:
                Bulk Velocity: list of doubles, length 3 (default: [0.0, 0.0, 0.0]) #[m/s]
                Temperature: double #[K]
                Number Density: double #[-/m^3]
```

**``Species``**: A list of the species that will have the given initial conditions.
Each entry needs to be present in the top-level ``Species`` block.

**``Constant``**: The parameters for a constant initial condion are specified under this key.
Incompatible with ``Sod``.

**``Bulk Velocity``**: The bulk velocity of the fluids in meters per second.

**``Temperature``**: The temperature in Kelvin of each fluid.
This is required to be positive.

**``Number Density``**: The number density of each fluid in number per meter cubed.
This is required to be positive.

**``Sod``**: The parameters for a Sod shock tube problem are specified under this key.
Incompatible with ``Constant``.

**``Discontinuity Location``**: The x location of the jump/discontinuity.

**``Left State``**: The state of the fluid when x is less than the discontinuity location.

**``Right State``**: The state of the fluid when x is greater than or equal to the discontinuity location.

#### Example

```yaml
Euler Fluids:
    Initial Conditions:
        - Species: [electron]
          Constant:
              Bulk Velocity: [1., 2., 4.]
              Temperature: 300
              Number Density: 1e20
        - Species: [electron, proton]
          Sod:
              Discontinuity Location: 0.7
              Left State:
                  Temperature: 20
                  Number Density: 1e22
              Right State:
                  Temperature: 320
                  Bulk Velocity: [-1., 0., 0.]
                  Number Density: 1e23
```

### Euler Fluids Boundary Conditions
The boundary conditions for the Euler fluids are specified under ```Boundary Conditions```.
The ```Boundary Conditions``` must be a sequence.
Each member of ```Boundary Conditions``` must specify the side for the boundary condition and the type of boundary condition to
apply.
The specified boundary conditions will be applied to all Euler fluids.

```yaml
Euler Fluids:
    Boundary Conditions:
        - Side: string or integer
          Type: string
```

**``Side``**: The side on which to apply the boundary condition.
For inline meshes the side should be a string.
The available sides are "left", "right", "top", "bottom", "front", and "back".
For meshes loaded from files the side should be a positive integer.
The integer is a mesh "attribute" in MFEM parlance or the index of the sideset in Cubit terminology.

**``Type``**: The type of particle boundary.
The options are:
- **``Reflecting``**: Reflects incident fluid.

Periodic boundaries are specified under ```Mesh``` and are handled automatically without any additional input here.

#### Example

```yaml
Euler Fluids:
    Boundary Conditions:
        - Side: left
          Type: Reflecting
        - Side: right
          Type: Reflecting
```

## Output

Output can be controlled under the ``Output`` key.
Two forms of output are produced:

- `Mesh Output`: this contains all fields defined on the mesh (e.g., electromagnetic fields and low-fidelity states). It comes in ParaView format.
- `Particle Dumps`: this contains the particle data (e.g., locations of particles in phase space). It comes in an hdf5 format.

The following options are available. All are optional.

### Stride
By specifying a ``Stride``, the user can control how often output is saved.
The ``Stride`` is the number of time steps between output dumps.
The default is `10`.

### Particle Dump Filename
The user can specify a preferred filename for the ``Particle Dumps``.
The default is ``particles.h5part``. If the ``h5part`` extension is not present, it is added.

### Mesh Output Folder
The user can specify a location to store the ``Mesh Output``.
The default is ``MeshOutput`` in the current working directory.

#### Example

```yaml
Output:
    Stride: 5
    Particle Dump Filename: my_particles.h5part
    Mesh Output Folder: MyMeshOutputs
```
