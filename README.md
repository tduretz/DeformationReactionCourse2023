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

The equilibrium condition can be reached iteratively. 

### Pseudo-transient integration

Here we use a pseudo-transient integration scheme and we state that:

$$
\frac{\partial H}{\partial \tau} = R
$$

The above equation can integrated in *pseudo time* ($\tau$) in an explicity manner such that:

$$
H^\mathrm{iter} = H^{\mathrm{iter}-1} + \Delta \tau R,
$$

where $\Delta \tau$ is the pseudo transient time step. It is a numerical parameter whose value is bounded by [stability analysis](https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis). $\mathrm{iter}$ corresponds to the pseudo transient iteration count. Once pseudo-transient steady state is achived, $\frac{\partial H}{\partial \tau} \rightarrow 0$, thus  $R \rightarrow 0$ and our target is reached: $\frac{\partial H}{\partial t} + \frac{\partial q}{\partial x} = 0$

### Pseudo-transient integration with damping



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

### 1D Simple shear model with fluid weakening

From the original [publication](https://www.sciencedirect.com/science/article/abs/pii/S0040195121003085) of Bras et al., 2021.

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_04.png" width=350px height=350px>

