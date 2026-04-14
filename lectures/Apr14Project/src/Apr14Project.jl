module Apr14Project

"""
    explicit_euler(f, u0::TR, t0, tmax, N) where TR <: Real

Advance the scalar initial value problem
``u'(t) = f(t, u(t))`` on ``[t0, tmax]`` using the explicit (forward) Euler method
with ``N`` uniform time steps.

Arguments:
- `f`: Right-hand side function with signature `f(t, u)`.
- `u0`: Scalar initial condition at time `t0`.
- `t0`: Initial time.
- `tmax`: Final time.
- `N`: Number of Euler steps (must be positive).

Returns:
- `(t_vals, u_vals)` where `t_vals` contains the time grid and `u_vals` the
    corresponding numerical solution values.

Notes:
- This is a first-order method with global error \$O(h)\$, where
    \$h = (tmax - t0)/N\$.
"""
function explicit_euler(f, u0::TR, t0, tmax, N) where TR <: Real
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

"""
    explicit_euler(f, u0::Vector{TR}, t0, tmax, N) where TR <: Real

Advance the vector-valued initial value problem
``u'(t) = f(t, u(t))`` on ``[t0, tmax]`` using the explicit (forward) Euler
method with ``N`` uniform time steps.

Arguments:
- `f`: Right-hand side function with signature `f(t, u)` returning a vector
    with the same dimension as `u`.
- `u0`: Vector initial condition at time `t0`.
- `t0`: Initial time.
- `tmax`: Final time.
- `N`: Number of Euler steps (must be positive).

Returns:
- `(t_vals, u_vals)` where `t_vals` is the time grid and `u_vals` is a vector
    of state vectors, one for each time point.

Notes:
- This is a first-order method with global error \$O(h)\$, where
    \$h = (tmax - t0)/N\$.
"""
function explicit_euler(f, u0::Vector{TR}, t0, tmax, N) where TR <: Real
    # step size
    h = (tmax - t0) / N

    # allocate arrays
    t_vals = Vector{Float64}(undef, N + 1)
    # preallocate a vector of vectors
    # u_vals = Vector{typeof(u0)}(undef, N + 1)
    u_vals = Vector{Vector{Float64}}(undef, N + 1)

    # initial condition
    t_vals[1] = t0
    u_vals[1] = copy(u0)

    # forward Euler iteration
    for n in 1:N
        t_vals[n+1] = t_vals[n] + h
        u_vals[n+1] = u_vals[n] + h * f(t_vals[n], u_vals[n])
    end

    return t_vals, u_vals
end


export explicit_euler
end # module Apr14Project
