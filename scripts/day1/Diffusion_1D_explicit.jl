using Plots

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
    nvx  = 100   # vertices
    ncx  = nvx-1 # centroids
    nt   = 1000  # time steps
    nout = 10    # plot each nout

    # Pre-compute 
    Δx   = (xmax-xmin)/ncx
    xv   = LinRange(xmin, xmax, nvx)
    xc   = 0.5*(xv[1:end-1] .+ xv[2:end])
    Δt   = Δx^2/(2.1*κ)

    # Allocate arrays
    T    = zeros(ncx)
    T0   = zeros(ncx+2) #--> diff(T0): ncx+2-1 = nvx
    qx   = zeros(nvx)   #--> diff(qx): nvx-1   = ncx
    Tana = zeros(ncx) # analytical solution

    # Initial condition
    T   .= Tmax.*exp.(-xc.^2/2/σ^2)
    Tini = copy(T)

    # Time loop
    for it=1:nt
        println( "step $it")
        # Populate T0
        T0[2:end-1] .= T
        T0[1]        = 2*Twest - T0[2]     # see section about Dirichlet BCs
        T0[end]      = 2*Teast - T0[end-1] # see section about Dirichlet BCs
        # Flux
        qx          .= -k.*diff(T0)/Δx
        # Explicit update
        T          .+= - Δt/(ρ*c) .* diff(qx)/Δx
        # Analytics
        t            = it*Δt
        Tana        .= Tmax ./ sqrt.(1 + 2*t*κ/σ^2) .* exp.(-xc.^2. / (2*(σ^2 + 2*t*κ)) )
        # Visualisation
        p = plot(xc,Tini, label="T initial")
        p = plot!(xc,T, label="T final")
        display(plot(p))
        sleep(0.01)
    end    
end

explicit_diffusion()