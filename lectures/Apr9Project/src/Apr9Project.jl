module Apr9Project


using Roots
# greet() = print("Hello World!")

function explicit_euler(f, u0, t0, tmax, N)
    # step size
    h = (tmax - t0) / N

    # allocate arrays
    t_vals = Vector{Float64}(undef, N + 1)
    u_vals = Vector{Float64}(undef, N + 1)

    # initial condition
    t_vals[1] = t0
    u_vals[1] = u0

    # forward Euler iteration
    for n in 1:N
        t_vals[n+1] = t_vals[n] + h
        u_vals[n+1] = u_vals[n] + h * f(t_vals[n], u_vals[n])
    end

    return t_vals, u_vals
end

function implicit_euler(f, u0, t0, tmax, N)
    h = (tmax - t0) / N

    # allocate arrays
    t_vals = Vector{Float64}(undef, N + 1)
    u_vals = Vector{Float64}(undef, N + 1)

    # initial condition
    t_vals[1] = t0
    u_vals[1] = u0

    for n in 1:N
        t_vals[n+1] = t_vals[n] + h

        # define the implicit Euler residual:
        # G(u) = u - u_n - h*f(t_{n+1}, u)
        G(u) = u - u_vals[n] - h * f(t_vals[n+1], u)

        # use previous value as the initial guess
        u_vals[n+1] = find_zero(G, u_vals[n])
    end

    return t_vals, u_vals
end

export explicit_euler, implicit_euler

end # module Apr9Project
