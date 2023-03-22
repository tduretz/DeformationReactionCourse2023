# Deformation Reaction Course 2023

The course takes place at Goethe University of Frankfurt. It is supervised by P. Yamato and T. Duretz, thanks to the support of the Heareus foundation.

|       |   |
| ----------- | ----------- |
| <img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/GU.svg" width=150px height=150px>       | <img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/Heraeus.png" width=150px>       |


## Day 0 - Get started

Once [Julia](https://julialang.org) is installed. From package mode (type ']'), you can activate the project (from the current directory) by typing:

`activate .`

If successful, you should see this: `(DeformationReactionCourse2023) pkg>`. In order to download/install all necessary packages required for this project, type:

`instantiate`

and wait a bit... If some packages fail to precompile, just restart Julia and then it should work.

## Day 1: Training

We get familiar with damped pseudo-transient integration for solving transient and steady state balance equations.

### Problem statement

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

where $\Delta \tau$ is the pseudo transient time step. It is a numerical parameter whose value is bounded by [stability analysis](https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis). $\mathrm{iter}$ corresponds to the pseudo transient iteration count. Once pseudo-transient steady state is achieved, $\frac{\partial H}{\partial \tau} \rightarrow 0$, thus  $R \rightarrow 0$ and our target is reached: $\frac{\partial H}{\partial t} + \frac{\partial q}{\partial x} = 0$

### Pseudo-transient integration with damping

The convergence of the above describe solution procedure can be drastically improved by accound for the *second order pseudo time derivative* of the quantitiy of interest:

$$
\theta_2 \frac{\partial^2 H}{\partial \tau^2} + \theta_1 \frac{\partial H}{\partial \tau} = R,
$$

where $\theta_1$ and $\theta_2$ are numerical parameters that control to the relative contribution of each derivative. This can be formulated as a 2 step solution method that includes the update of the pseudo rate:

$$
\frac{\partial H}{\partial \tau}^\mathrm{iter}= R + (1-\theta_1)\frac{\partial H}{\partial \tau}^{\mathrm{iter}-1}
$$

and the subsequent update of the quantity:

$$
H^\mathrm{iter} = H^{\mathrm{iter}-1}  + \Delta \tau \frac{\partial H}{\partial \tau}^\mathrm{iter}.
$$

Here we have ome more crucial numerical parameter $\theta_1$. Efficient convergence can be achieved by fine tuning this parameter. There is so far no automatic way to determine the optimal value of this parameter. Analytical values can be obtained for idealised cases are were proven to be relatively robust for some cases ([Raess et al., 2022](https://gmd.copernicus.org/articles/15/5757/2022/)). 
To our knowledge, the earliest description of this solution method can be found in [Frankel, 1950](https://www.jstor.org/stable/2002770). One main advantage of this approach is the possibility to port codes to GPU computing in a rather simple manner, very useful resources canbe found on the [PDE on GPU website](https://pde-on-gpu.vaw.ethz.ch). 

### 1D heat equation 

We seek to integrate the following equations ($T$):

$$
\begin{align}
q = -k\frac{\partial{T}}{\partial{x}} \\
\rho c\frac{\partial{T}}{\partial{t}}=-\frac{\partial{q}}{\partial{x}} \\
\end{align}
$$

Flowchart of the pseudo-transient iteration cycle:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_02.png" width=350px>

Each application of a `diff` operator reduces array size by `-1`.

## Day 2 - 1D Simple shear model

In this part, we solve the following equation for velocity ($v_x$).

$$
\begin{align}
\tau_{xy} = \eta \frac{\partial v_x}{\partial y} \\
0=\frac{\partial{\tau_{xy}}}{\partial y}
\end{align}
$$

Variable arrangement for the 1D mechanical problem:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_01.png" width=350px>

Flowchart of the pseudo-transient iteration cycle:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_03.png" width=350px>

Each application of a `diff` operator reduces array size by `-1`.

## Day 2 - 1D Simple shear model with fluid weakening

We now couple a 1D mechanics with the diffusion equation for fluid amount:

$$
\begin{align}
\eta &= f(F) \\
\tau_{xy} &= \eta \frac{\partial v_x}{\partial y} \\
q_F &= -D  \frac{\partial F}{\partial y}    \\
0 &= \frac{\partial{\tau_{xy}}}{\partial y} \\
\frac{\partial{F}}{\partial{t}} &= - \frac{\partial q_F}{\partial y}
\end{align}
$$

This model assumes that the material is ready to be transformed (i.e., pressure and temperature are sufficiently large) and that only fluid is missing to start transformation. Once fluid diffuses through the material, it initiates transformation and viscosity drops, which triggers the development of a shear zone. This model is an crude approximation of a multi-phase model, thus the evolution of porosity and fluid pressure is not modelled. This is taken from the original [publication](https://www.sciencedirect.com/science/article/abs/pii/S0040195121003085) of Bras et al., 2021.

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_04.png" width=350px>

More implementation details

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_05.png" width=350px height=350px>

Example of solution:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/ShearZoneFluidWeakening.gif" width=350px>

