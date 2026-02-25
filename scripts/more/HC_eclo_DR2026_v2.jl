# 1D HC model - Eclogitisation of granulite rock (Bras et al., 2023) 
using Plots, LinearAlgebra, StaticArrays, ForwardDiff
using SpecialFunctions
using Printf
using MAT

const derivative = true
const primitive  = false 

function update_Pf(Pf, Pf_phys, Pf_it, ρT, t, p, LU)

    rel = p.rel

    ρT[1]   = t < p.t_pert ?  2*p.ρTWest-ρT[2] : ρT[2]
    ρT[end] = ρT[end-1]

    for i=1:size(Pf,1)-0
        Pf_phys[i] = Itp1D_rev_scalar1(LU.ρT_lt, LU.Pf_lt, ρT[i])
        Pf[i]      = rel*Pf_phys[i] + (1-rel)*Pf_it[i]
    end

end

function residual_local(ρT, ρT0, Pf_it, xc, t, p, LU, Δx, Δt, i, ncx)
    rel = p.rel

    # BC on ρT ?
    if i==2
        ρT[1] = t < p.t_pert ?  2*p.ρTWest-ρT[2] : ρT[2]
    end
    if i==ncx-1
        ρT[3] = ρT[2]
    end

    # Database
    Pf_phys = SVector{3}( Itp1D_rev_scalar1(LU.ρT_lt, LU.Pf_lt, ρT[ii]) for ii in 1:3)
    Pf  = SVector{3}( rel*Pf_phys[ii]+(1-rel)*Pf_it[ii]  for ii in 1:3) 
    ρs  = SVector{3}( Itp1D_scalar1( LU.Pf_lt, LU.ρs_lt[:], Pf[ii], LU.ΔP_lt, LU.Pf_lt[1]) for ii in 1:3)
    ρf  = SVector{3}( Itp1D_scalar1( LU.Pf_lt, LU.ρf_lt[:], Pf[ii], LU.ΔP_lt, LU.Pf_lt[1]) for ii in 1:3)

    Xs  = p.Xs_g  .- (p.Xs_g .- p.Xs_e) .* (Pf .> p.Pf_gr_out)
    
    ϕ   = 1.0 .- p.ρs_tot ./ ρs ./ Xs 

    # Averaging
    ϕv  = @. 0.5 * (ϕ[1:end-1]  + ϕ[2:end])
    ρfv = @. 0.5 * (ρf[1:end-1] + ρf[2:end])

    # Darcy flux:
    k_ηf     = @. p.k_ηf0 * ϕv^p.npow

    # Total mass flux
    qρT      = @. -ρfv * k_ηf * (Pf[2:end] - Pf[1:end-1])/Δx

    # Residual
    Rρ       = -(ρT[2]  - ρT0 ) / Δt - (qρT[2] - qρT[1]) / Δx  
    # if t < p.t_pert && xc <= p.wini
    #     Rρ       = -(ρT[2] ) / Δt /10000 
    # end

    return Rρ
end

function residual(r, r0, ρT, ρT0, Pf_it, xc, t, p, LU, Δx, Δt, deriv)
    for i=2:size(r,1)-1 
        ρT_loc    = @MVector [ρT[i-1], ρT[i], ρT[i+1]]
        Pf_it_loc = @MVector [Pf_it[i-1], Pf_it[i], Pf_it[i+1]]
        if deriv == false
            r0[i] = r[i]
            r[i]  = residual_local(ρT_loc, ρT0[i], Pf_it_loc, xc[i], t, p, LU, Δx, Δt, i, length(xc))
        else
            ∂r∂P = ForwardDiff.gradient(x->residual_local(x, ρT0[i], Pf_it_loc, xc[i], t, p, LU, Δx, Δt, i, length(xc)), ρT_loc) 
            r0[i] = abs(∂r∂P[2])
            r[i]  = sum(abs.(∂r∂P))
        end
    end
end

@views function Itp1D_rev_scalar1(xlt, varlt, xdata)
    xinf_id = sum( xlt .- xdata .< 0 )
    if xinf_id<1 xinf_id=1 end 
    if xinf_id>length(xlt)-1 xinf_id=length(xlt)-1 end 
    xinf_dist = (xlt[xinf_id+1] - xdata) / (xlt[xinf_id+1] - xlt[xinf_id]);
    return xinf_dist * varlt[xinf_id] + (1. - xinf_dist) * varlt[xinf_id+1];
end

@views function Itp1D_scalar1(xlt, varlt, xdata, dx, xmin)
    iW = Int(floor((xdata - xmin) / dx) + 1)
    if iW<1 iW=1 end
    if iW>size(xlt,1)-2 iW=size(xlt,1)-2 end
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
    Xs_g      = 1.0         # Prop. of the solid phase in granulite
    Xs_e      = 0.993       # Prop. of the solid phase in eclogite

    # Model numerics : --------------------------------------
    CFL_t    = 0.09
    rel      = 1e-2
    eps_err  = 1e-5
    ncx      = 101                    # number of grid points
    nout     = 50
    nt       = 1000
    dt_fact  = 100
    niter    = 1e5
    # preprocessing
    Δx       = Lx/ncx              # grid spacing
    xc       = LinRange(-Δx/2, Lx+Δx/2, ncx+2)  # grid points cooarrdinates
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
    Pf      = zeros(ncx+2)
    Pf_it   = zeros(ncx+2)
    Pf_phys = zeros(ncx+2)
    ρs      = zeros(ncx+2)
    ρf      = zeros(ncx+2)
    ρT      = zeros(ncx+2)
    ρT0     = zeros(ncx+2)
    βf      = zeros(ncx+2)
    Xs      = zeros(ncx+2)
    ϕ       = zeros(ncx+2)
    dρT_dτ  = zeros(ncx+2)
    Rρ      = zeros(ncx+2)
    # DR
    D       =  ones(ncx+2)
    G       =  ones(ncx+2)
    Rρ0     = zeros(ncx+2)
    # Allocate arrays for vertices
    ϕv      = zeros(ncx+1)
    k_ηf    = zeros(ncx+1)
    # Initialise
    Pf     .= Pbg        # Background pressure without fluid pressure perturbation
    Pf[xc .<= wini] .+= Pamp
    Pfi     = copy(Pf)
    PfWest  = Pfi[1]
    # Density look up - Initial perturbation is fully eclogitised from the start
    for ip = 1:ncx+2
        ρs[ip] = Itp1D_scalar1( Pf_lt, ρs_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
        ρf[ip] = Itp1D_scalar1( Pf_lt, ρf_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
        βf[ip] = Itp1D_scalar1( Pf_lt, βf_lt[:], Pf[ip], ΔP_lt, Pf_lt[1])
    end
    # Xs = (1.0 - XH2O)
    Xs  .= Xs_g  .- (Xs_g .- Xs_e) .* (Pf .> Pf_gr_out)
    # Porosity:
    ϕ   .= 1.0 .- ρs_tot ./ ρs ./ Xs
    # Total density:
    ρT  .= (1.0 .- ϕ) .* ρs .+ ϕ .* ρf

    p = (
        Pf_gr_out = Pf_gr_out,
        PfWest    = PfWest,
        ρTWest    = ρT[1],
        npow      = npow,
        k_ηf0     = k_ηf0,
        wini      = wini,
        t_pert    = t_pert,
        Xs_g      = Xs_g,
        Xs_e      = Xs_e,
        ρs_tot    = ρs_tot,
        rel       = rel,
    )

    LU = (
        ρT_lt = ρT_lt, 
        Pf_lt = Pf_lt,
        ΔP_lt = ΔP_lt,
        ρs_lt = ρs_lt,
        ρf_lt = ρf_lt,
        βf_lt = βf_lt,
    )

    Pf_it .= Pf

    # TIME LOOP:
    anim = @animate for it = 1:nt

        # Old values
        ρT0 .= ρT

        # Define time step
        ϕv   .= 0.5* (ϕ[1:end-1] + ϕ[2:end])
        k_ηf .= k_ηf0 .*ϕv .^npow
        Dcmax = maximum(max(k_ηf[1:end-1],k_ηf[2:end]) ./ βf[2:end-1])
        Δτ0   = CFL_t*Δx^2*minimum(ρf)/maximum(ρT)/Dcmax
        Δt    = Δτ0*dt_fact

        # Define PT params
        CFL    = 0.9
        c_fact = 0.9
        nrρ0   = 1.0
        residual(G, D, ρT, ρT0, Pf_it, xc, t, p, LU, Δx, Δt, derivative)
        λmax   = maximum(G ./ D)
        Δτ     = 2 / sqrt(maximum(λmax)) * CFL
        λmin   = 0.0
        c      = 2*sqrt(λmin)*c_fact
        α      = 2 * Δτ^2 / (2 + c.*Δτ)
        β      = (2 - c * Δτ) / (2 + c.*Δτ)

        # Pseudo transient loop
        for iter = 1:niter

            Pf_it .= Pf

            # Residual
            residual(Rρ, Rρ0, ρT, ρT0, Pf_it, xc, t, p, LU, Δx, Δt, primitive)

            # Update pseudo-rate
            dρT_dτ .=  Rρ ./ D .+ β .* dρT_dτ

            # Update density
            ρT    .+=  α * dρT_dτ

            # Update pressure
            update_Pf(Pf, Pf_phys, Pf_it, ρT, t, p, LU)

            if iter==1 || mod(iter, nout) == 0

                # Exit ?
                nrρ     = norm(Rρ)/length(Rρ)
                if iter==1 nrρ0 = nrρ; end
                @printf("Step %04d --- Iteration %05d\n", it, iter)
                @printf("||rρ|| = %2.10e --- ||rρ/rρ0|| = %2.10e\n", nrρ, nrρ/nrρ0)
                @show norm(Pf.-Pf_phys)

                if (min(nrρ, nrρ/nrρ0) < eps_err) break; end

                # Define PT params
                residual(G, D, ρT, ρT0, Pf_it, xc, t, p, LU, Δx, Δt, derivative)
                λmax   = maximum(G ./ D)
                Δτ     = 2 / sqrt(maximum(λmax)) * CFL    
                λmin   = abs.((sum(Δτ.*dρT_dτ.*( (Rρ .- Rρ0) ./ D )))) / sum( (Δτ.*dρT_dτ).^2 )
                c      = 2 * sqrt(λmin) * c_fact
                α      = 2 * Δτ^2 / (2 + c*Δτ)
                β      = (2 - c * Δτ) / (2 + c*Δτ)
            end
        end

        # Final updates
        t = t + Δt

        # Plot results:
        if it==1 || mod(it,5) == 0
            tD = t/24/3600
            p1 = plot(xc, Rρ)
            p2 = plot(xc, Pfi./1e9,color=:black, label= "Pfi", marker = 0.0,linewidth = 1.0,
                    xlabel = "xc [m]", ylabel = "P [GPa]")
            p2 = plot!(xc, 1.94.*ones(size(xc)), linestyle=:dash, color=:green, label= "Pr", marker = 0.0,linewidth = 2.0,α=0.5)
            p2 = plot!(xc, Pf./1e9,color=:red, linestyle=:dash, label= "Pf", 
                    marker = 2.0,linewidth = 1.0,
                    title = " t = $(round(tD,digits=2)) days")
            display(plot(p2))
        end
    end
end

main_HC()