# 1D shear-zone widening associated with eclogitization reaction (Bras et al., 2021) 
using Plots, SpecialFunctions, Printf
import LinearAlgebra:norm

function main_shear_zone_fluid()
    # Physical parameters:
    ymin       = 0.0       # ymin [m]
    ymax       = 2.0       # ymax [m]
    ε̇BG        = 1e-14     # background strainrate [s-1]
    Wini       = 0.05      # width of initially hydrated zone
    D          = 1e-14     # Fluid diffusivity [m2/s]
    Fini       = 100.0     # initial amount of fluid  
    Vnorth     = 2.0 * ε̇BG * (ymax-ymin)
    Vsouth     = 0.0
    # Material viscosity 
    ηi         = 1.0e22    # strong initial phase
    ηf         = 1.0e19    # final phase
    # Numerical parameters:
    nvy        = 100       # nb of nodes
    ncy        = nvy-1     # nb of cells
    nt         = 100        # nb of physical time steps
    niter      = 10000      # nb of iterations
    nout       = 100       # outputs in iterations 
    rel_tol_V  = 1e-5      # tolerance for velocity 
    rel_tol_F  = 1e-5      # tolerance for water diffusion
    # --------------------
    Ly         = ymax-ymin
    Δy         = Ly/ncy
    yv         = LinRange(ymin, ymax, nvy)
    yce        = LinRange(ymin-Δy/2, ymax+Δy/2, ncy+2)    
    # For velocity
    Vx         = zeros(ncy+2)
    Ux         = zeros(nvy+1) # displacement
    Ry         = zeros(ncy)
    dVxdτ      = zeros(ncy)  # for damping
    ηv         = zeros(nvy)
    τxy        = zeros(nvy)
    # Fluid
    F          = zeros(ncy+2)
    F0         = zeros(ncy+2)
    qF         = zeros(nvy)
    Rf         = zeros(ncy)
    dFdτ       = zeros(ncy)
    # Reaction
    r          = zeros(ncy+2)
    ϕ          = zeros(ncy+2) 
    ϕv         = zeros(nvy)
    # Initial conditions
    t               = 0.
    F[yce .< Wini] .= Fini
    ηv             .= ηi
    r[yce .< Wini] .= 1.0
    # Numerics
    Δt         = Δy^2/D/2.1 
    Δτv        = Δy^2/4.1/maximum(ηv)  # pseudo-transient time step for velocity
    Δτx        = 2*Δt*Δy^2/(4.1*D*Δt + Δy^2)
    θv         = Ly/nvy * 0.5
    θf         = Ly/nvy * 4.3
    # Physical time loop
    for it = 1:nt
        # Initialise time step
        t  += Δt
        F0 .= F 
        Vx .= 0.0
        rv0, rf0 = 1.0, 1.0
        # Pseudo transient iterations
        for iter = 1:niter
            # (0) Set boundary conditions
            F[1]         = F[2]                  # zero flux
            F[end]       = F[end-1]
            Vx[1]        = 2*Vsouth - Vx[2]      # Dirichlet
            Vx[end]      = 2*Vnorth - Vx[end-1] 
            # Non linearity
            ϕv          .= 0.5.*(ϕ[1:end-1] .+ ϕ[2:end])
            ηv          .= ((1.0 .- ϕv)./ηi .+ ϕv./ηf).^(-1)
            Δτv = Δy^2/4.1/maximum(ηv)  # change PT time step
            # (1) strain rate, stress for stokes & flux for X
            qF  .= -D .* diff(F, dims=1)  ./ Δy
            τxy .= ηv .* diff(Vx, dims=1) ./ Δy
            # (3) residual 
            Ry .=    diff(τxy, dims=1) ./ Δy
            Rf .= .- diff(qF , dims=1) ./ Δy .- (F[2:end-1] .- F0[2:end-1]) ./ Δt
            # (4) update pseudo transient rates
            dVxdτ .= Ry .+ (1. - θv) .* dVxdτ  # with damping
            dFdτ  .= Rf .+ (1. - θf) .* dFdτ   # with damping
            # (5) update Vx and X
            Vx[2:end-1] .+= Δτv .* dVxdτ
            F[2:end-1]  .+= Δτx .* dFdτ
            # Relate viscosity to ϕ
            ϕ .= F
            ϕ[F .>= 1.0] .= 1.0
            # Checks
            if iter==1 || mod(iter,nout) == 0
                rv     = norm(Ry)/length(Ry)
                rf     = norm(Rf)/length(Rf)
                if iter==1 rf0 = rf; rv0 = rv; end
                @printf("it %04d --- iteration %05d --- ||rv|| = %2.10e\n", it, iter, rv/rv0)
                @printf("it %04d --- iteration %05d --- ||rf|| = %2.10e\n", it, iter, rf/rf0)
                if (rv/rv0<rel_tol_V && rf/rf0<rel_tol_F) break; end
            end
        end 
        # Update displacement
        Ux .+= Vx .* Δt 
        # Figures
        plot()
        p1 = plot(Ux, yce, color=:red,label="" ,marker = 0,linewidth = 1.0,
            xlabel = "displacement [m]",
            ylabel = "y [m]",xlims=(0.0, 0.2),ylims=(0,2.0))
        p2 = plot(ϕ, yce, color=:blue,label="" ,marker = 0,linewidth = 2.0,
            xlabel = "ϕ",
            ylabel = "y [m]", xlims=(0,1), ylims=(0,2.0))
        p3 = plot(r, yce, color=:black, label="" , marker = 0, linewidth = 2.0,
            xlabel = "reaction",
            ylabel = "y [m]", xlims=(0,1.01),ylims=(0,2.0))
        p4 = plot(log10.(ηv), yv, color=:black,label="", marker = 0, linewidth = 2.0,
            xlabel = "log10(η) Pa.s",
            ylabel = "y [m]", xlims=(18.9,22.1), ylims=(0,2.0))
        p = plot(p1, p2, p3, p4, layout = (2,2))          
        display(p)
        sleep(0.01)
    end
end

main_shear_zone_fluid()