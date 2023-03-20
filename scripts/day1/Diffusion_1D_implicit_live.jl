using Plots, Printf, LinearAlgebra

function explicit_diffusion()

    # Physics
    xmin  = -1/2
    xmax  = 1/2
    Tmax  = 100   # max. amplitude
    σ     = 0.1   # bandwiΔth
    k     = 1.    # thermal conductivity
    ρ     = 1000  # density
    c     = 1000  # heat capacity
    κ     = k/ρ/c # diffusivity
    Twest = 0.0 
    Teast = 0.0  

    # Numerics
    nvx   = 100    # vertices
    ncx   = nvx-1  # centroids
    nt    = 100    # time steps
    niter = 100000 # number of iterations
    nout  = 10     # plot each nout

    # Pre-compute 
    Δx   = (xmax-xmin)/ncx
    xv   = LinRange(xmin, xmax, nvx)
    xc   = 0.5*(xv[1:end-1] .+ xv[2:end])
    Δt   = Δx^2/(2.1*κ) 
    Δτ   = 2*Δt*Δx^2/(4.1*κ*Δt + Δx^2) 

    # Allocate arrays
    T    = zeros(ncx+2) 
    R    = zeros(ncx)
    T0   = zeros(ncx+2) # diff(T0) --> ncx+2-1 = nvx
    qx   = zeros(nvx)   # diff(qx) --> nvx-1   = ncx
    Tana = zeros(ncx)   # analytical solution

    # Initial condition
    T[2:end-1] .= Tmax.*exp.(-xc.^2/2/σ^2)
    Tini        = copy(T)

    # Time loop
    for it=1:nt
        # Get solution from previous time step
        T0 .= T
        # Show step number
        println( "step $it")
        # Pseudo-time loop (iterations)
        for iter=1:niter
            # Populate T
            T[1]        = 2*Twest - T0[2]     # see section about Dirichlet BCs
            T[end]      = 2*Teast - T0[end-1] # see section about Dirichlet BCs
            # Flux
            qx          .= -k.*diff(T)/Δx
            # Residual
            R           .= ρ*c*(T[2:end-1] .- T0[2:end-1])./Δt .+ diff(qx)./Δx
            # Explicit update
            T[2:end-1] .-= Δτ/(ρ*c) .* R
            # Monitor error
            err          = norm(R)/length(R)
            @printf("it = %04d ---- iter = %05d --- err = %2.4e\n", it, iter, err)
            if err<1e-6
                break
            end
        end
        # Analytics
        t            = it*Δt
        Tana        .= Tmax ./ sqrt.(1 + 2*t*κ/σ^2) .* exp.(-xc.^2. / (2*(σ^2 + 2*t*κ)) )
        # Visualisation
        if mod(it, nout)==0 || it==1
            p = plot( xc, Tini[2:end-1], label="T initial")
            p = plot!(xc, T[2:end-1], label="T final")
            p = plot!(xc, Tana, label="T anal")
            display(plot(p))
            sleep(0.1)
        end
    end    
end

explicit_diffusion()