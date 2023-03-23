## Day 3

Check out [pdf presentation](https://hessenbox-a10.rz.uni-frankfurt.de/getlink/fi9gsHsvSbfsihdGAtW44d/ReactionDef_Part02_compressed.pdf).
This is taken from the original publication of [Yamato et al., 2022](https://www.sciencedirect.com/science/article/abs/pii/S0012821X2200156X).

The governing equations are:

$$
\begin{align}
\nabla{v} = \frac{1}{r} \frac{\partial(r v_r)}{\partial r} \\
\dot{\varepsilon}_{rr}=\frac{\partial{v_r}}{\partial r} - \frac{1}{3}\nabla{v}\\
\dot{\varepsilon}_{\phi\phi}=\frac{v_{r}}{r} - \frac{1}{3}\nabla{v}\\
\dot{\varepsilon}_{rr}= \frac{\tau_{rr}}{2\eta} + \frac{1}{2G}\frac{\mathrm{d} \tau_{rr}}{\mathrm{d}  t}\\
\dot{\varepsilon}_{\phi \phi}= \frac{\tau_{\phi \phi}}{2\eta} + \frac{1}{2G}\frac{\mathrm{d}  \tau_{\phi \phi}}{\mathrm{d} t}\\
\sigma_{rr}=-P + \tau_{rr} \\
\sigma_{\phi \phi}=-P + \tau_{\phi \phi} \\
0=\frac{\partial \sigma_{rr}}{\partial r}+\frac{1}{r}\left( \sigma_{rr}-\sigma_{\phi\phi}\right) \\
\frac{\mathrm{d} \ln (\rho)}{\mathrm{d}  t} = -\nabla v\\
\end{align}
$$

The density is made function of reaction progress ($X$) depending on a kinetic law and presure ($P$):

$$
\begin{align}
ρ = ρ_0 \exp (β P)\\
ρ_0 = Xρ_0^{f} + (1−X)ρ_0^{i} \\
\frac{\mathrm{d} X}{\mathrm{d} t}=\frac{X_{eq}-X}{\tau_k}\\
X_{eq}=1-\frac{1}{2} \left( \mathrm{erfc}\left( \frac{P-P_r}{d_{P_r}}\right)\right)
\end{align}
$$

Computation of $X$ using kinetic law:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_06.png" width=350px>

Computation of $\rho_0$ from $X$ and computation of $\rho$ from $\rho_0$ and $P$:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/pict_07.png" width=350px>

Evolution of density using the kinetic law:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/X_density_evolution.gif" width=350px>

Summary of relations that need to be evaluated:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/step01_Equations.png" width=350px>

Implementation of Maxwwell viscoselastic model:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/step03_Elasco_viscosity.png" width=350px>

Application to $rr$ and $\phi\phi$ deviatoric stress components:

<img src="https://github.com/tduretz/DeformationReactionCourse2023/blob/main/images/step04_Final_Maxwell_Formulation.png" width=350px>
