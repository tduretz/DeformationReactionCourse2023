using Plots, Printf
import LinearAlgebra: norm

function main_shear()

# Physics
ymin   = -5.0
ymax   = 5.0
η      = 1e20
Vnorth = 1.0
Vsouth = -1.0  

# Numerics
ny   = 100  # vertices
ncy  = ny-1 # centroids
nout = 10   # plot each nout

# Pre-compute
Δy   = (ymax-ymin)/ncy 
yv   = LinRange(ymin, ymax, ny)
yce  = LinRange(ymin-Δx/2, ymax+Δy/2, ncy+2)
Δτ   = Δy^2/2/η/2.1
θ    = 0.6*(ymax-ymin)/ncy

# Allocate arrays
Vx    = zeros(ncy+2)
τxy   = zeros(ny)
dVxdτ = zeros(ncy)
R     = zeros(ncy)
err0  = 0.
# Pseudo-time loop
for iter=1:10000
    # Set boundary conditions
    Vx[1]         = 2*Vsouth - Vx[2]
    Vx[end]       = 2*Vnorth - Vx[end-1]
    # Flux
    τxy          .= 2*η.*diff(Vx, dims=1)/Δy
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

plot()
plot!(yce, Vx, label="Vx", xlabel="Vx", ylabel="y")

end

main_shear()