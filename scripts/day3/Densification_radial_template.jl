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

function compDensity!(ρ, ρ0, X, Xo, Xeq, P, Δt, Pr, dPr, τk, ρ0i, ρ0f, β)
    Xeq .= 1 .- 0.5 .* (erfc.((P .- Pr)./dPr));
    X   .= Δt/(τk+Δt) .* Xeq .+ τk/(τk+Δt) .* Xo;
    ρ0  .= X .* ρ0f .+ (1 .-X) .* ρ0i
    ρ   .= ρ0 .* exp.(β .* P)
end

function main()
    # Spatial domain:
    rmin    = 0.0        # m
    rmax    = 1.1284     # m
    rinc    = 0.25       # m
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
    nt        = 100        # nb of steps
    Δt        = τk/10      # t step
    ηe        = G*Δt       # ηe definition
    ncr       = 100        # number of cells 
    niter     = 1e5        # nb of iterations
    # Pseudo-transient parameters
    αV        = 0.03 
    αP        = 3.00
    θ         = 5*(rmax-rmin)/ncr
    nout      = 5000
    rel_tol_V = 1e-6
    rel_tol_P = 1e-6
    # Define grid with variable spacing
    rc, rv, Δrc, Δri, rv0, ncr_inc = InitialiseGrid(ncr, rmax, rmin, rinc)
    # ...    
end # End function

main()
