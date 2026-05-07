module May7Project

using SparseArrays

function backwards_difference_matrix(nx, dx)
    nx >= 1 || throw(ArgumentError("nx must be at least 1"))
    dx != 0 || throw(ArgumentError("dx must be nonzero"))

    # explictly constructs an nx by nx matrix of floating point numbers
    Dx = zeros(nx, nx)
    # builds the matrix row by
    Dx[1, 1] = 1 / dx
    for j in 2:nx
        Dx[j, j] = 1 / dx
        Dx[j, j - 1] = -1 / dx
    end

    Dx
end

function sparse_backwards_difference_matrix(nx, dx)
    nx >= 1 || throw(ArgumentError("nx must be at least 1"))
    dx != 0 || throw(ArgumentError("dx must be nonzero"))

    diag = fill(1 / dx, nx)
    subdiag = fill(-1 / dx, nx - 1)

    spdiagm(-1 => subdiag, 0 => diag)
end

function sparse_backwards_difference_matrix_pbc(nx, dx)
    nx >= 1 || throw(ArgumentError("nx must be at least 1"))
    dx != 0 || throw(ArgumentError("dx must be nonzero"))

    diag = fill(1 / dx, nx)
    subdiag = fill(-1 / dx, nx - 1)

    Dx = spdiagm(-1 => subdiag, 0 => diag)
    Dx[1, nx] = -1 / dx

    Dx
end

# code to be compatible with DifferentialEquations.jl
function advection_mf!(du, u, p, t)
    # unpack paramters
    c = p[1]; # c in the PDE
    dx = p[2]; # mesh spacing
    x = p[3]; # mesh points where we are solving
    g = p[4]; # boundary condition function at x = a
    f = p[5]; # source function in the PDE

    du[1] = -c/dx * (u[1] - g(t)) + f(x[1], t)
    for j in 2:length(u)
        du[j] = -c/dx * (u[j] - u[j-1]) + f(x[j], t)
    end
    
    du

end

function advection!(du, u, p, t)
    # unpack paramters
    c = p[1]; # c in the PDE
    dx = p[2]; # mesh spacing
    x = p[3]; # mesh points where we are solving
    g = p[4]; # boundary condition function at x = a
    f = p[5]; # source function in the PDE
    Dx = p[6]; # differentiation matrix for the spatial discretization

    # .= assigns the values into the exisitn array du
    du .= -c * Dx * u + f.(x, t)
    du[1] += c/dx *g(t); # enforce the boundary condition at x = a

    
    du

end


function sparse_second_derivative_matrix(n, h)
    @assert n >= 1 "n must be at least 1"
    @assert h > 0 "h must be positive"

    main = fill(-2.0, n)
    off  = fill(1.0, n - 1)

    return (1 / h^2) * spdiagm(-1 => off, 0 => main, 1 => off)
end

export backwards_difference_matrix, sparse_backwards_difference_matrix, sparse_backwards_difference_matrix_pbc, advection_mf!, advection!, sparse_second_derivative_matrix


end # module May7Project
