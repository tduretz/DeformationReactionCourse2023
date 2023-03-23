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
