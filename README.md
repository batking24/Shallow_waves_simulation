# Shallow Water Equations Simulation

This repository contains a Python implementation of the shallow water equations (SWE) simulation. The code simulates fluid dynamics in a shallow water environment, providing both 2D vector field visualization and 3D wave animation of the fluid behavior.

## Introduction

The shallow water equations are a set of hyperbolic partial differential equations that describe the flow of a fluid under the effects of gravity, where the horizontal length scale is much greater than the vertical length scale. These equations are derived from the Navier-Stokes equations, which govern fluid motion, by depth-averaging under the assumption that vertical acceleration is negligible.

### Physical Background

The simulation is based on three main equations:
1. Conservation of mass (continuity equation)
2. Conservation of momentum in x-direction
3. Conservation of momentum in y-direction

In our implementation, we use:

```
du/dt - fv = -g*d(eta)/dx + tau_x/(rho_0*H) - kappa*u
dv/dt + fu = -g*d(eta)/dy + tau_y/(rho_0*H) - kappa*v
d(eta)/dt + d((eta + H)*u)/dx + d((eta + H)*u)/dy = sigma - w
```

Where:
- u, v: velocity components in x and y directions
- η (eta): surface elevation
- f: Coriolis parameter
- g: gravitational acceleration
- H: mean water depth
- τ: wind stress
- ρ₀: reference density
- κ: bottom friction coefficient

## Features

- Simulation of shallow water dynamics with:
  - Coriolis force
  - Bottom friction
  - Wind stress
  - Variable bathymetry
- 2D vector field visualization of fluid velocities
- 3D animation of surface elevation
- Customizable physical and computational parameters
- Time performance analysis
- CFL (Courant-Friedrichs-Lewy) condition implementation for numerical stability

## Dependencies

- NumPy: For numerical computations
- Matplotlib: For visualization and animation
- mpl_toolkits.mplot3d: For 3D surface plots

## Code Structure

The code is organized into several main sections:

1. **Visualization Tools**
   - `velocity_animation`: Creates 2D vector field animations
   - `eta_animation3D`: Generates 3D surface elevation animations

2. **Parameter Settings**
   - Physical parameters (e.g., gravity, depth, Coriolis parameter)
   - Computational parameters (grid size, time steps)
   - Boolean flags for different physical effects

3. **Array Initialization**
   - Setup of grid points and meshgrid
   - Initialization of velocity and elevation arrays
   - Temporary arrays for upwind scheme calculations

4. **Main Computation Loop**
   - Implementation of the shallow water equations
   - Upwind scheme for numerical stability
   - Boundary conditions
   - Time stepping

5. **Visualization**
   - Vector field plotting
   - 3D surface animation
   - Data sampling for animation

## Key Parameters

```python
# Physical Parameters
L_x = 1E+6          # Domain length in x-direction [m]
L_y = 1E+6          # Domain length in y-direction [m]
g = 9.81            # Acceleration of gravity [m/s²]
H = 100             # Depth of fluid [m]
f_0 = 1E-4         # Fixed part of Coriolis parameter [1/s]
beta = 2E-11       # Gradient of Coriolis parameter [1/ms]
rho_0 = 1024.0     # Density of fluid [kg/m³]

# Computational Parameters
N_x = 150          # Number of grid points in x-direction
N_y = 150          # Number of grid points in y-direction
max_time_step = 5000  # Total number of time steps
```

## Numerical Methods

### Upwind Scheme
The code implements an upwind scheme for solving the hyperbolic PDEs, which is crucial for numerical stability. This method uses the direction of fluid flow to determine the spatial discretization.

### CFL Condition
The time step is chosen according to the Courant-Friedrichs-Lewy (CFL) condition:
```python
dt = 0.1*min(dx, dy)/np.sqrt(g*H)
```
This ensures numerical stability by making sure that information doesn't propagate faster than the numerical scheme can handle.

## Output Examples

### 2D Vector Field Visualization
<video src="./velocity_compress.mp4" controls="controls" style="max-width: 730px;">
</video>

This visualization shows the velocity field of the fluid, with arrows indicating direction and magnitude of flow.

### 3D Wave Animation
<video src="./eta_surface_compressed.mp4" controls="controls" style="max-width: 730px;">
</video>




The 3D animation displays the evolution of surface elevation over time, showing wave propagation and interaction.

## Performance Analysis

For default parameters (N_x = N_y = 150):
- Visualization tools: 74 ms
- Parameter setting: 56 ms
- Array initialization: 68 ms
- Main computation: 11.1 s
- Vector field visualization: 40 ms
- 3D wave visualization: 2min 53sec

For increased resolution (N_x = N_y = 200):
- Visualization tools: 70 ms
- Parameter setting: 56 ms
- Array initialization: 95 ms
- Main computation: 17.1 s
- Vector field visualization: 1min 53sec
- 3D wave visualization: 2min 53sec

## Usage

1. Ensure all dependencies are installed
2. Adjust parameters as needed in the code
3. Run the simulation
4. View the generated animations

## Future Improvements

Potential areas for enhancement:
- Parallel computation for improved performance
- Additional physical effects (e.g., temperature, salinity)
- Interactive parameter adjustment
- Real-time visualization
- Implementation of different boundary conditions

## License

[Insert your chosen license information here]

## Contributors

[Your name/organization]

## References

1. Pedlosky, J. (1987). Geophysical Fluid Dynamics. Springer
2. Cushman-Roisin, B. (1994). Introduction to Geophysical Fluid Dynamics
3. Numerical Recipes in Scientific Computing
