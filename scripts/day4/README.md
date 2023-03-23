The aim of the exercise if to couple reaction and fluid flow. For this purpose, data have been preliminary generated using equilibrium themodynamic computations (e.g. [PerpleX](https://www.perplex.ethz.ch)). They are referred to as look-up tables ('_lt'). The main partial differential equation is the total mass conservation, which can be solved for the total density ($\rho_T$). Other variables are read (i.e, *interpolated*) from the look-up table (`LUT_plagio_eclo.mat`). The governing equations are formulated as:

$$
\begin{align}
\rho_f=f(P) \\
\rho_s=f(P) \\
\beta_f=f(P) \\
P_f = f(\rho_T) \\
\phi=1-\frac{\rho_T}{\rho_s X_s} \\
q_{\rho_T} = -\rho_f\frac{k\phi^3}{\eta_f}\left( \frac{\partial P_f}{\partial x}\right)\\
\frac{\partial \rho_T}{\partial t}=-\frac{\partial q}{\partial x}
\end{align}
$$

[Link to the lecture](https://hessenbox-a10.rz.uni-frankfurt.de/getlink/fiKVUb5ZAUMuQsTjgBEHUA/ReactionDef_Part03_compressed.pdf)

Computation of solid fraction $X_s$, porosity $\phi$ and total density $\rho^T$:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_08.png" width=500px>

Porosity only evolves due to the transition between granulte and eclogite:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_09.png" width=500px>

Discretisation:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_10.png" width=500px>
