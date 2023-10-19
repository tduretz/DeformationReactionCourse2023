## Day 1: Training

We get familiar with damped pseudo-transient integration for solving transient and steady state balance equations. [Lecture notes](https://hessenbox-a10.rz.uni-frankfurt.de/getlink/fiXUcFqSEWx22HVLHjom1z/2_SolvingEquations_compressed.pdf).

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

Here we have one more crucial numerical parameter $\theta_1$. Efficient convergence can be achieved by fine tuning this parameter. There is so far no automatic way to determine the optimal value of this parameter. Analytical values can be obtained for idealised cases are were proven to be relatively robust for some cases ([Raess et al., 2022](https://gmd.copernicus.org/articles/15/5757/2022/)). 
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
