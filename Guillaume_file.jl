using Plots, Printf, LinearAlgebra
using SpecialFunctions

# Generate the main function
function explicit_diffusion()

    # Physical parameters
    xmin = -1 / 2
    xmax = 1 / 2
    Tmax = 100 # Initial amplitude
    σ = 0.1 # bandwidth
    k = 1.0 # thermal conductivity
    ρ = 1000 # density
    cp = 1000 # heat capacity
    κ = k / (ρ * cp)
    Twest = 0.0
    Teast = 0.0

    # Numerical paramaters
    nvx = 100      # vertices
    ncx = nvx - 1  # centroids
    nt = 1  # time steps
    niter = 100
    nout = 10  # plot each out

    # Pre-compute
    Δx = (xmax - xmin) / ncx
    xv = LinRange(xmin, xmax, nvx)
    xc = 0.5 * (xv[1:end-1] .+ xv[2:end])
    Δt = Δx^2 / (κ * 2.1)
    Δτ = 2 * Δt * Δx^2 / (4 * κ * Δt + Δx^2)

    # Allocate space to store values of T
    T = zeros(ncx + 2)
    R = zeros(ncx)
    T0 = zeros(ncx + 2)
    qx = zeros(nvx)
    Tana = zeros(ncx)

    # Initial condition
    T[2:end-1] .= Tmax .* exp.(-xc .^ 2 / 2 / σ^2)
    Tini = copy(T)

    # Time loop
    for it = 1:nt

        # Get solution from previous time step
        T0 .= T

        # Show step number
        # println("step $it")

        # Pseudo time loop (iterations)
        for iter = 1:niter

            # Populate T0
            T[1] = 2 * Twest - T0[2]
            T[end] = 2 * Teast - T0[end-1]

            # Flux
            qx .= -k .* diff(T0, dims=1) / Δx

            # Residuals
            R .= ρ * cp * (T[2:end-1] .- T0[2:end-1]) ./ Δt .+ diff(qx) ./ Δx

            # Explicit update
            T[2:end-1] .-= Δτ / (ρ * cp) .* R

            # Monitor error
            err = norm(R) / length(R)
            @printf("it = %04d -- iter = %05d -- error = %2.4e\n", it, iter, err)
            if err < 1e-6
                break
            end
        end

        # Analytics
        t = it * Δt
        Tana .= Tmax ./ sqrt.(1 + 2 * t * κ / σ^2) .* exp.(-xc .^ 2 / (2 * (σ^2 + 2 * t * κ)))

        # p = plot(xc, Tini, label="T initial")
        # p = plot!(xc, T[2:end-1], label="T $it")
        # p = plot!(xc, Tana[2:end-1], label="T analytical $it")
        # display(plot(p))
        # sleep(0.001)
    end

end # End main function

explicit_diffusion()
