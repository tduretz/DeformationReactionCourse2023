# Shrinking 1D (Yamato et al., 2022 - EPSL)
using Plots
using SpecialFunctions
using Printf
import LinearAlgebra:norm
kyr     = 1000*365.25*24*3600

function InitialiseGrid(ncr, rmax, rmin, rinc)
    Lr      = rmax - rmin   # Model length
    L_inc   = rinc - rmin
    L_mat   = rmax - rinc
    ncr_inc = Int(floor(L_inc/Lr*ncr))
    ncr_mat = ncr - ncr_inc
    Δr_inc  = L_inc/ncr_inc
    Δr_mat  = L_mat/ncr_mat
    Δrc     = [ones(ncr_inc).* Δr_inc; ones(ncr_mat).* Δr_mat]
    Δri     = 0.5 * (Δrc[2:end] + Δrc[1:end-1])
    rv, rv0 = [0; cumsum(Δrc)] , [0; cumsum(Δrc)]
    rc      = 0.5 * (rv[2:end] + rv[1:end-1])
    return rc, rv, Δrc, Δri, rv0, ncr_inc
end

function main()
    # Spatial domain:
    rmin      = 0.0        # m
    rmax      = 1.1284     # m
    rinc      = 0.25       # m
    # material properties (see Table 1):
    ρ0i       = 2850       # kg.m-3
    ρ0f       = 3250       # kg.m-3
    Pr        = 1.5e9      # Pa
    dPr       = 700e6      # Pa
    Pbg       = 2.0e9      # Pa
    K         = 80e9       # Pa
    β         = 1/K        # Pa-1  
    τk        = 3.1558e9   # s
    ηv        = 1e22       # Pa.s
    G         = 40e9       # Pa
    # Numerical parameters
    nt        = 1       # nb of steps
    Δt        = τk/10      # t step
    ηe        = G*Δt       # ηe definition
    nvr       = 101        # number of cells
    ncr       = nvr-1      # number of cells 
    niter     = 1e5        # nb of iterations
    # Pseudo-transient parameters
    αV        = 0.03 
    αP        = 3.00
    θ         = 1.4*(rmax-rmin)/ncr
    nout      = 5000
    rel_tol_V = 1e-6
    rel_tol_P = 1e-6
    # Define grid with variable spacing
    rc, rv, Δrc, Δri, rv0, ncr_inc = InitialiseGrid(ncr, rmax, rmin, rinc)
    # Allocate arrays
    Rv, ∂V∂τ      = zeros(nvr), zeros(nvr)
    Rp, ∂P∂τ      = zeros(ncr), zeros(ncr)
    Vr, Vr0, Vrc  = zeros(nvr), zeros(nvr), zeros(ncr)
    ∇V            = zeros(ncr)
    P             = zeros(ncr)
    ρ, ρ0,  ρr    = zeros(ncr), zeros(ncr), zeros(ncr)
    τrr, τrr0     = zeros(ncr), zeros(ncr)
    τϕϕ, τϕϕ0     = zeros(ncr), zeros(ncr)
    σrr           = zeros(ncr)
    σrri          = zeros(size(rc,1)-1)
    σϕϕ           = zeros(ncr)
    σϕϕi          = zeros(size(rc,1)-1)
    τzz           = zeros(ncr)
    τII           = zeros(ncr)
    ε̇rr           = zeros(ncr)
    ε̇ϕϕ           = zeros(ncr)
    Ėrr           = zeros(ncr)
    Ėϕϕ           = zeros(ncr)
    phc           = zeros(ncr)
    X, X0         = zeros(ncr), zeros(ncr)
    # Initialisation
    Vr[2:end-1]     .= 1e-20.*rand(nvr-2)
    P               .= Pbg
    ρ               .= ρ0i .* exp.(β*P)
    phc[rc .> rinc] .= 1    
    t                = 0.0
    # Storage for plots:
    r_inc_vec  = zeros(nt)
    p_inc_vec  = zeros(nt)
    p_mat_vec  = zeros(nt)
    t_vec      = zeros(nt)
    maxτII_vec = zeros(nt)
    X_inc_vec  = zeros(nt)
    # TIME LOOP
    for it = 1:nt
        # Initialise t step
        t    += Δt
        # Store old value
        τrr0 .= τrr
        τϕϕ0 .= τϕϕ
        ρ0   .= ρ
        X0   .= X
        # For error monitoring
        rr0, rp0 = 0., 0.
        # Pseudo transient iterations
        for iter = 1:niter
            # Kinematics:
            Vrc   .= 0.5 * (Vr[1:end-1] + Vr[2:end])
            ∇V    .= 1.0 ./rc .* diff(rv.*Vr, dims=1) ./ Δrc      
            ε̇rr   .= diff(Vr,dims=1) ./ Δrc - 1/3 .* ∇V
            ε̇ϕϕ   .= Vrc./rc  - 1/3 .* ∇V
            Ėϕϕ   .= ε̇ϕϕ .+ τϕϕ0./2.0./ηe
            Ėrr   .= ε̇rr .+ τrr0./2.0./ηe    
            # Rheology:
            ηve    = (1.0/ηe + 1.0/ηv)^(-1)
            τrr   .= 2.0*ηve*Ėrr
            τϕϕ   .= 2.0*ηve*Ėϕϕ
            σrr   .= -P .+ τrr
            σϕϕ   .= -P .+ τϕϕ
            # Non-linearity
            X     .= (X0 .* τk + (1.0 .- 0.5 .* erfc.((P .- Pr) ./ dPr)) .* Δt) ./ (τk + Δt)
            X[phc.==1] .= 0 
            ρr    .= X.*ρ0f .+ (1.0 .-X).*ρ0i
            ρ     .= ρr .*exp.( β .*P)
            # Interpolation
            σrri  .= 0.5 .* (σrr[1:end-1] .+ σrr[2:end])
            σϕϕi  .= 0.5 .* (σϕϕ[1:end-1] .+ σϕϕ[2:end])
            # Residuals
            Rv[2:end-1] .=  diff(σrr,dims=1)./Δri .+ 1.0 ./ rv[2:end-1] .* (σrri.-σϕϕi)
            Rp          .= .- ( (log.(ρ) .- log.(ρ0))./Δt .+ ∇V )
            # Checks
            if iter==1 || mod(iter,nout) == 0
                # Ėrrors
                rr     = norm(Rv)/length(Rv)
                rp     = norm(Rp) /length(Rp)
                if iter==1 rr0 = rr; rp0 = rp; end
                @printf("Iteration %05d\n", iter)
                @printf("||fr|| = %2.10e\n", rr/rr0)
                @printf("||fp|| = %2.10e\n", rp/rp0)
                if (rr/rr0<rel_tol_V && rp/rp0<rel_tol_P) break; end
            end
            # PT steps
            ΔtV    = αV*minimum(Δrc)^2/ηve/2.1
            ΔtP    = αP*ηve/ncr/2.1
            # Rate updates
            ∂V∂τ   .= Rv .+ (1.0-θ) .* ∂V∂τ 
            ∂P∂τ   .= Rp
            # Solution updates
            Vr    .+= ΔtV.*∂V∂τ
            P     .+= ΔtP.*∂P∂τ        
        end
        if it == 1
            rv0 .= rv
            Vr0 .= Vr
        end
        # Update mesh
        rv  .= rv .+ Δt .* Vr
        rc  .= 0.5*(rv[2:end] + rv[1:end-1])
        Δrc .= rv[2:end] - rv[1:end-1]
        Δri .= rc[2:end] - rc[1:end-1]
        # Compute deviators
        τzz .= -τϕϕ.-τrr
        τII .= sqrt.(0.5 .* ( τrr.^2 .+ τzz.^2 .+ τϕϕ.^2))
        # Store solutions
        r_inc_vec[it]  = rv[ncr_inc+1]
        p_inc_vec[it]  = P[1]
        p_mat_vec[it]  = P[end]
        t_vec[it]   = t
        maxτII_vec[it] = maximum(τII)
        X_inc_vec[it]  = X[1]
        # Visualisation
        if mod(it, 10)==0
            plot()
            pvecplot = LinRange(0.0,3.0,100)
            Xvecplot = zeros(size(pvecplot))
            Xvecplot .= 1.0 .- 0.5.*erfc.((pvecplot.*1e9 .- Pr)./dPr)
            p1 = plot!(pvecplot, Xvecplot, marker = 0, 
                    xlabel = "P [GPa]",
                    ylabel = "X")
            # Subplot 1
            p1 = plot!(p_inc_vec[1:it]./1e9,X_inc_vec[1:it], color=:red, label="", marker=0, linewidth=1.0)
            p1 = plot!(p_mat_vec[1:it]./1e9,X_inc_vec[1:it], color=:blue,label="", marker=0, linewidth=1.0)
            # Subplot 2
            p2 = plot( t_vec[1:it]./kyr, p_inc_vec[1:it]./1e9, color=:green, label="", marker=0, linewidth=1.0)
            p2 = plot!(t_vec[1:it]./kyr, p_mat_vec[1:it]./1e9, color=:black, label="", marker=0, linewidth=1.0,
                        xlabel = "t [kyrs]",
                        ylabel = "P [GPa]")
            # Subplot 3
            p3 = plot(t_vec[1:it]./kyr, maxτII_vec[1:it]./1e9, color=:red,label="" , marker=0, linewidth=1.0,
                        xlabel = "t [kyrs]",
                        ylabel = "max, τII [GPa]")
            p = plot(p1,p2,p3, layout = (3,1))   
            display(p)
        end
    end
end # End function

main()