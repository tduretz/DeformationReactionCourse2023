# 1D HC model - Eclogitisation of granulite rock (Bras et al., 2023) 
using Plots
using SpecialFunctions
using Printf
using MAT

@views function Itp1D_rev_scalar1(xlt, varlt, xdata)
    xinf_id = sum( xlt .- xdata .< 0 )
    if xinf_id<1 xinf_id=1 end 
    if xinf_id>length(xlt)-1 xinf_id=length(xlt)-1 end 
    xinf_dist = (xlt[xinf_id+1] - xdata) / (xlt[xinf_id+1] - xlt[xinf_id]);
    return xinf_dist * varlt[xinf_id] + (1. - xinf_dist) * varlt[xinf_id+1];
end

@views function Itp1D_scalar1(xlt, varlt, xdata, dx, xmin)
    iW = Int(floor((xdata - xmin) / dx) + 1)
    wW = 1. - (xdata - xlt[iW]) / dx
    return wW * varlt[iW] + (1. - wW) * varlt[iW + 1]
end

function LoadData()
    # Read Lookup tables: ----------------------------------- 
    file      = matopen( string(@__DIR__, "/LUT_plagio_eclo.mat") )
    Pf_lt     = Array(read(file, "Pf_lt"))[:]
    ρg_lt     = Array(read(file, "rho_g"))[:]
    ρe_lt     = Array(read(file, "rho_e"))[:]
    ρw_lt     = Array(read(file, "rho_w"))[:]
    n_lu      = 2001
    close(file) 
    ρs_lt = zeros(n_lu)
    ρf_lt = zeros(n_lu)
    ρT_lt = zeros(n_lu)
    Xs_lt = zeros(n_lu)
    ϕ_lt  = zeros(n_lu)
    # Infos from LT:
    ΔP_lt      = Pf_lt[2] - Pf_lt[1]
    Pfc_lt     = 0.5* (Pf_lt[1:end-1] + Pf_lt[2:end])/2
    ieclo_in   = findfirst(xc-> xc > 0.0 ,ρe_lt)  # first nodes of Pf_lt in eclogite
    Pf_eclo_in = Pf_lt[ieclo_in]                # corresponding Pf value
    Pf_gr_out  = Pf_lt[ieclo_in-1]              # previous Pf value
    # DENSITIES and COMPRESSIBILITY:
    # ρsolid = ρgranulite when Pf< Pr and ρeclogite when Pf>Pr
    ρs_lt .= ρg_lt                                
    ρs_lt[Pf_lt .>= Pf_eclo_in] .= ρe_lt[Pf_lt .>= Pf_eclo_in]
    # ρfluid = ρwater 
    ρf_lt .= ρw_lt
    # get β:
    βf_lt = diff(log.(ρf_lt), dims=1) ./ diff(Pf_lt, dims=1)
    return Pf_lt, Pfc_lt, ΔP_lt, βf_lt, ρf_lt, ρs_lt, ρT_lt, Xs_lt, ϕ_lt, ρg_lt, ρe_lt, ρw_lt, Pf_gr_out
end

function main_HC()

    # Physical parameters: -----------------------------------
    Lx        = 0.05        # model length [m]
    Pbg       = 15.0e8      # 
    Pamp      = 1.0e9       # Pressure perturbation (1 GPa)
    wini      = 0.005       # Initial perturbation width (half of initial width)
    t_pert    = 0.3*24*3600 # t of the fluid pulse
    k_ηf0     = 1e-19/1e-3  # m^2/s/Pa - k = 1e-19 m2, etaf = 1e3 Pa.s
    npow      = 3;          #
    Xs_g      = 1.0              # Prop. of the solid phase in granulite
    Xs_e      = 0.993            # Prop. of the solid phase in eclogite

    # Model numerics : --------------------------------------
    CFL      = 0.09
    rel0     = 0.25
    eps_err  = 1e-3
    ncx      = 101                    # number of grid points
    nout     = 10
    nt       = 1000
    θ        = 0.75  
    dt_fact  = 100
    niter    = 1e5
    # preprocessing
    Δx       = Lx/ncx              # grid spacing
    xc       = LinRange(0.0, Lx, ncx)  # grid points cooarrdinates
    t        = 0.0
    Δt       = 0.0
    # Read Lookup tables: ----------------------------------- 
    Pf_lt, Pfc_lt, ΔP_lt, βf_lt, ρf_lt, ρs_lt, ρT_lt, Xs_lt, ϕ_lt, ρg_lt, ρe_lt, ρw_lt, Pf_gr_out = LoadData() 

    # Total solid mass (ρgranulite @ initial pressure Pbg)---
    # ---> Interpolate ρs_tot = f(P_inf) from the database 
    ρs_tot = Itp1D_scalar1( Pf_lt, ρg_lt, Pbg, ΔP_lt, Pf_lt[1])

    # proportion of the solid phase:
    # ---> Compute Xs_lt using Xs_g, Xs_e, Xs_g
    # ---> Compute ϕ_lt using Xs_lt
    # ---> Compute ρT_lt using ϕ_lt, ρf_lt
    Xs_lt  .= Xs_g .* ones(size(Pf_lt)) .- (Xs_g .- Xs_e) .* (Pf_lt .> Pf_gr_out)
    ϕ_lt   .= 1.0 .- ρs_tot ./ ρs_lt ./ Xs_lt
    ρT_lt  .= (1.0 .- ϕ_lt) .* ρs_lt .+ ϕ_lt .* ρf_lt
    # Visualisation
    # Plot ρs_tot = f(Pf_lt)
    # Plot ρs_lt  = f(Pf_lt)
    # Plot ρf_lt  = f(Pf_lt)
    # Plot ρT_lt  = f(Pf_lt)
    plot()
    p1 = plot!(Pf_lt./1e9, ρs_tot .* ones(size(Pf_lt)),color=:orange, linestyle=:dash, label= "ρStot", marker = 0.0,linewidth = 2.0)
    p1 = plot!(Pf_lt./1e9, ρs_lt,color=:cyan, label= "ρs(LUT)", marker = 0.0,linewidth = 2.0)
    p1 = plot!(Pf_lt./1e9, ρf_lt,color=:red, label= "ρf(LUT)", marker = 0.0,linewidth = 2.0)
    p1 = plot!(Pf_lt./1e9, ρT_lt,color=:black, linestyle=:dash, label= "ρT(LUT)", marker = 0.0,linewidth = 2.0)
    display(p1)
    # Allocate centroid arrays
    Pf     = zeros(ncx)
    Pfr    = zeros(ncx)
    ρs     = zeros(ncx)
    ρf     = zeros(ncx)
    ρT     = zeros(ncx)
    ρT0    = zeros(ncx)
    βf     = zeros(ncx)
    Xs     = zeros(ncx)
    ϕ      = zeros(ncx)
    dρT_dt = zeros(ncx)
    dρT_dτ = zeros(ncx)
    Rρ     = zeros(ncx)
    # Allocate extended centroids arrays
    ϕe     = zeros(ncx+2)
    ρfe    = zeros(ncx+2)
    Pfe    = zeros(ncx+2)
    # Allocate arrays for vertices
    ϕv     = zeros(ncx+1)
    ρfv    = zeros(ncx+1)
    k_ηf   = zeros(ncx+1)
    qρT    = zeros(ncx+1)
    # Initialise
    Pf   .= Pbg        # Background pressure without fluid pressure perturbation
    Pf[xc .<= wini] .+= Pamp
    Pfi    = copy(Pf)
    PfWest = Pfi[1]
    # Density look up - Initial perturbation is fully eclogitised from the start
    for ip = 1:ncx
        ρs[ip] = Itp1D_scalar1( Pf_lt, ρs_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
        ρf[ip] = Itp1D_scalar1( Pf_lt, ρf_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
        βf[ip] = Itp1D_scalar1( Pf_lt, βf_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
    end
    # Xs = (1.0 - XH2O)
    Xs  .= Xs_g .* ones(ncx) .- (Xs_g .- Xs_e) .* (Pf .> Pf_gr_out)
    # Porosity:
    ϕ   .= 1.0 .- ρs_tot ./ ρs ./ Xs
    # Total density:
    ρT  .= (1.0 .- ϕ) .* ρs .+ ϕ .* ρf

    # TIME LOOP:
    anim = @animate for it = 1:nt
        # Old values
        rel  = rel0
        ρT0 .= ρT
        # Relative residual
        rρ0  = 1.0
        # Define transient time step: ------
        ϕe   .= [ϕ[1]; ϕ[1:end]; ϕ[end]]
        ϕv   .= 0.5* (ϕe[1:end-1] + ϕe[2:end])
        k_ηf .= k_ηf0 .*ϕv .^npow
        Dcmax = maximum(max(k_ηf[1:end-1],k_ηf[2:end]) ./ βf)
        Δτ    = CFL*Δx^2*minimum(ρf)/maximum(ρT)/Dcmax
        Δt    = Δτ*dt_fact
        # Pseudo transient loop
        for iter = 1:niter
            # (0) Reverse LUT : Get fluid pressure:
            for ip=1:ncx
                Pfr[ip] = Itp1D_rev_scalar1(ρT_lt, Pf_lt, ρT[ip])
            end
            Pf .= (1.0 .- rel) .* Pf .+ rel .* Pfr
            # (1) Density look up + Xs [already done before] ---
            for ip = 1:ncx 
                ρs[ip] = Itp1D_scalar1( Pf_lt, ρs_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
                ρf[ip] = Itp1D_scalar1( Pf_lt, ρf_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
                βf[ip] = Itp1D_scalar1( Pf_lt, βf_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
            end
            Xs  .= Xs_g .* ones(ncx) .- (Xs_g .- Xs_e) .* (Pf .> Pf_gr_out)
        
            # (2) Get porosity as fct of ρs and Xh20 (equ. 13) ---
            ϕ     .= 1.0 .- ρs_tot ./ ρs ./ Xs        
            
            # boundary conditions and averaging
            ϕe  .= [ϕ[1]; ϕ; ϕ[end]]
            ρfe .= [ρf[1]; ρf; ρf[end]]

            if t < t_pert
                 Pfe .= [2.0 .* PfWest .- Pf[2]; Pf[1:end]; Pfe[end-1]]
            else
                 Pfe .= [Pf[2]; Pf; Pf[end-1]]
            end
            
            ϕv  .= 0.5* (ϕe[1:end-1] + ϕe[2:end])
            ρfv .= 0.5* (ρfe[1:end-1] + ρfe[2:end])

            # (3) total mass conservation with darcy flux:
            k_ηf     .= k_ηf0 .*ϕv .^npow
            qρT      .= .-ρfv .* k_ηf .* diff(Pfe)./Δx

            # Residual
            Rρ      .= .-(ρT .- ρT0) ./ Δt .- diff(qρT)./Δx  
            if t < t_pert Rρ[xc .<= wini] .= 0.0; end
          
            # Update rate         
            dρT_dτ   .= Rρ .+ (1 - θ) .* dρT_dτ

            # Update rho
            ρT      .+= Δτ .* dρT_dτ

            if iter==1 || mod(iter,nout) == 0
                rρ     = norm(Rρ)/length(Rρ)
                if iter==1 rρ0=rρ; end
                @printf("Step %04d --- Iteration %05d\n", it, iter)
                @printf("||rρ|| = %2.10e\n", rρ/rρ0)
                if (rρ/rρ0 < eps_err) break; end
            end

        end

        t = t + Δt

        # Plot results:
        if it==1 || mod(it,5) == 0
            tD = t/24/3600
            p2 = plot(xc,Pfi./1e9,color=:black, label= "Pfi", marker = 0.0,linewidth = 1.0,
                    xlabel = "xc [m]", ylabel = "P [GPa]")
            p2 = plot!(xc,1.94.*ones(size(xc)), linestyle=:dash, color=:green, label= "Pr", marker = 0.0,linewidth = 2.0,α=0.5)
            p2 = plot!(xc,Pf./1e9,color=:red, linestyle=:dash, label= "Pf", 
                    marker = 2.0,linewidth = 1.0,
                    title = " t = $(round(tD,digits=2)) days")
            display(p2)
        end
    end
end

main_HC()