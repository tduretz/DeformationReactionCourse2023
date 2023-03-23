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
    ieclo_in   = findfirst(x-> x > 0.0 ,ρe_lt)  # first nodes of Pf_lt in eclogite
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
    Pbg       = 15.0e8           # Background pressure 
    Xs_g      = 1.0              # Prop. of the solid phase in granulite
    Xs_e      = 0.993            # Prop. of the solid phase in eclogite
    
    # Read Lookup tables: ----------------------------------- 
    Pf_lt, Pfc_lt, ΔP_lt, βf_lt, ρf_lt, ρs_lt, ρT_lt, Xs_lt, ϕ_lt, ρg_lt, ρe_lt, ρw_lt, Pf_gr_out = LoadData() 

    # Total solid mass (ρgranulite @ initial pressure Pbg)---
    # ---> Interpolate ρs_tot = f(P_inf) from the database 

    # proportion of the solid phase:
    # ---> Compute Xs_lt using Xs_g, Xs_e, Xs_g
    # ---> Compute ϕ_lt using Xs_lt
    # ---> Compute ρT_lt using ϕ_lt, ρf_lt


    # Visualisation
    # ---> Plot ρs_tot = f(Pf_lt)
    # ---> Plot ρs_lt  = f(Pf_lt)
    # ---> Plot ρf_lt  = f(Pf_lt)
    # ---> Plot ρT_lt  = f(Pf_lt)

end

main_HC()