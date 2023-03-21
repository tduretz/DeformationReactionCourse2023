using Plots, Printf
import LinearAlgebra: norm

function main_shear_variable_viscosity()

    # Physics
    ymin   = -0.5
    ymax   = 0.5
    η0     = 1e20 
    σ      = 0.1
    Vnorth = 1.0
    Vsouth = -1.0  

    # Numerics
    nvy    = 100  # vertices
    ncy    = nvy-1 # centroids
    nout   = 10   # plot each nout

    # Pre-compute
    Δy     = (ymax-ymin)/ncy 
    yv     = LinRange(ymin, ymax, nvy)
    yce    = LinRange(ymin-Δy/2, ymax+Δy/2, ncy+2)

    # Allocate arrays
    Vx     = zeros(ncy+2)
    τxy    = zeros(nvy)
    dVxdτ  = zeros(ncy)
    R      = zeros(ncy)

    # Viscosity 
    ηv     = zeros(nvy)
    ηv     = η0 .* (1  .- 0.9.* exp.(-yv.^2/2/σ^2))
    err0   = 0.

    Δτ     = Δy^2/2/maximum(ηv)/2.1
    θ      = 6*(ymax-ymin)/ncy

    # Pseudo-time loop
    for iter=1:10000
        # Set boundary conditions
        Vx[1]     = 2*Vsouth - Vx[2]
        Vx[end]   = 2*Vnorth - Vx[end-1]
        # Flux
        τxy      .= 2*ηv.*diff(Vx, dims=1)/Δy
        # Residual
        R        .= diff(τxy, dims=1)/Δy
        dVxdτ    .= R .+ (1.0 - θ).*dVxdτ
        # Check error
        err = norm(R)/length(R)
        if iter==1 err0 = err; end
        if (err/err0) < 1e-6 break; end
        @printf("iter %04d - error = %2.6e\n", iter, err/err0)
        # Update
        Vx[2:end-1]  .+= dVxdτ * Δτ
    end

    p1 = plot(Vx, yce, label="Vx final", xlabel="Vx", ylabel="y")
    p2 = plot(diff(Vx,dims=1)./Δy, yv, label="ε̇xy", xlabel="Vx", ylabel="y")
    display(plot(p1,p2))

end

main_shear_variable_viscosity()