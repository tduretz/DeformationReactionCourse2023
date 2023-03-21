# Deformation Reaction Course 2023

## Day 0 - Get started

Once [Julia](https://julialang.org) is installed. From package mode (type ']'), you can activate the project (from the current directory) by typing:

`activate .`

If successful, you should see this: `(DeformationReactionCourse2023) pkg>`. In order to download/install all necessary packages required for this project, type:

`instantiate`

and wait a bit... If some packages fail to precompile, just restart Julia and then it should work.

## Day 1 - Day 2: Training

We get familiar with damped pseudo-transient integration for solving transient and steady state balance equations.
The aim is to solve balance laws of the following form:

$$
\frac{\partial H}{\partial t} + \frac{\partial q}{\partial x} = 0
$$

where $H$ is the quantity of interest and $q$ is the its flux (e.g. advection, diffusion...). We seek for an implicit solution, the equilibrium can not be obtained directly, instead there will be an imbalance or residual, $R$:

$$
R = \frac{\partial H}{\partial t} + \frac{\partial q}{\partial x}.
$$

The equilibrium condition can be reached iteratively. Here we use a pseudo-transient integration scheme and we state that:

$$
\frac{\partial H}{\partial \tau} = R
$$

The above equation can integrated in *pseudo time* ($\tau$) in an explicity manner such that:

$$
H^\mathrm{new} = H^\mathrm{old} + \Delta \tau R,
$$

where $\Delta \tau$ is the pseudo transient time step. It is a numerical parameter whose value is bounded by [stability analysis](https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis).



### 1D heat equation 

Flowchart of the pseudo-transient iteration cycle:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_02.png" width=350px height=350px>

Each application of a `diff` operator reduces array size by `-1`.

### 1D Simple shear model

Variable arrangement for the 1D mechanical problem:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_01.png" width=350px height=350px>

Flowchart of the pseudo-transient iteration cycle:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_03.png" width=350px height=350px>

Each application of a `diff` operator reduces array size by `-1`.
