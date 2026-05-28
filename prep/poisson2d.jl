using SparseArrays

function assemble_laplacian2d(Δx, Δy, Nx, Ny)
    # Construct the 2D Laplacian for the rectangle with Dirichlet boundary conditions

    rows = Int[] # rows
    cols = Int[] # columuns
    vals = Float64[]  # values

    for k in 1:Nx*Ny
        # get local coordinates from global coordinate
        j, i = fldmod1(k, Nx)
        # corners, then edges, then interior
        if (i == 1 && j == 1)
            append!(rows, k * ones(Int, 3))
            append!(cols, [1, 2, 1 + Nx])
            append!(vals, [-2 / Δx^2 - 2 / Δy^2, 1 / Δx^2, 1 / Δy^2])
        elseif (i == 1 && j == Ny)
            append!(rows, k * ones(Int, 3))
            append!(cols, [1 + (Ny - 2) * Nx, 1 + (Ny - 1) * Nx, 2 + (Ny - 1) * Nx])
            append!(vals, [1 / Δy^2, -2 / Δx^2 - 2 / Δy^2, 1 / Δx^2])
        elseif (i == Nx && j == 1)
            append!(rows, k * ones(Int, 3))
            append!(cols, [Nx - 1, Nx, Nx + Nx])
            append!(vals, [1 / Δx^2, -2 / Δx^2 - 2 / Δy^2, 1 / Δy^2])
        elseif (i == Nx && j == Ny)
            append!(rows, k * ones(Int, 3))
            append!(cols, [Nx - 1 + Nx * (Ny-1), Nx + Nx * (Ny-2), Nx * Ny])
            append!(vals, [1 / Δx^2, 1 / Δy^2, -2 / Δx^2 - 2 / Δy^2])
        elseif (i == 1 && j > 1 && j < Ny)
            # right edge of domain, (1<j<Ny)
            append!(rows, k * ones(Int, 4))
            append!(cols, [1 + (j - 2) * Nx, 1 + (j - 1) * Nx, 2 + (j - 1) * Nx, 1 + j * Nx])
            append!(vals, [1 / Δy^2, -2 / Δx^2 - 2 / Δy^2, 1 / Δx^2, 1 / Δy^2])
        elseif (i == Nx && j > 1 && j < Ny)
            # left edge of domain, (1<j<Ny)
            append!(rows, k * ones(Int, 4))
            append!(cols, [(j - 1) * Nx, j * Nx - 1, j * Nx, (j + 1) * Nx])
            append!(vals, [1 / Δy^2, 1 / Δx^2, -2 / Δx^2 - 2 / Δy^2, 1 / Δy^2])
        elseif (i > 1 && i < Nx && j == 1)
            # bottom edge of domain, (1<i<Nx)
            append!(rows, k * ones(Int, 4))
            append!(cols, [i - 1, i, i + 1, i + Nx])
            append!(vals, [1 / Δx^2, -2 / Δx^2 - 2 / Δy^2, 1 / Δx^2, 1 / Δy^2])
        elseif (i > 1 && i < Nx && j == Ny)
            # top edge of domain, (1<i<Nx)
            append!(rows, k * ones(Int, 4))
            append!(cols, [i + (Ny - 2) * Nx, i - 1 + (Ny - 1) * Nx, i + (Ny - 1) * Nx, i + 1 + (Ny - 1) * Nx])
            append!(vals, [1 / Δy^2, 1 / Δx^2, -2 / Δx^2 - 2 / Δy^2, 1 / Δx^2])
        else
            # interior of domain
            append!(rows, (i + (j - 1) * Nx) * ones(Int, 5))
            append!(cols, [i + (j - 2) * Nx, i - 1 + (j - 1) * Nx, i + (j - 1) * Nx, i + 1 + (j - 1) * Nx, i + j * Nx])
            append!(vals, [1 / Δy^2, 1 / Δx^2, -2 / Δx^2 - 2 / Δy^2, 1 / Δx^2, 1 / Δy^2])
        end
    end

    A = sparse(rows, cols, vals, Nx * Ny, Nx * Ny)
    #     display(cols)
    return A
end


function dirichlet_laplacian_2d_kron(Δx, Δy, Nx, Ny)

    # number of interior points
    # Nx = nx - 1
    # Ny = ny - 1

    # 1D Dirichlet Laplacian in x
    ex = ones(Nx)
    Tx = spdiagm(
        -1 => ex[1:end-1],
         0 => -2 .* ex,
         1 => ex[1:end-1]
    ) / Δx^2

    # 1D Dirichlet Laplacian in y
    ey = ones(Ny)
    Ty = spdiagm(
        -1 => ey[1:end-1],
         0 => -2 .* ey,
         1 => ey[1:end-1]
    ) / Δy^2

    Ix = sparse(I, Nx, Nx)
    Iy = sparse(I, Ny, Ny)

    # 2D Laplacian
    L = kron(Iy, Tx) + kron(Ty, Ix)

    return L
end



function sparse_second_derivative_matrix(n, h)
    @assert n >= 1 "n must be at least 1"
    @assert h > 0 "h must be positive"

    main = fill(-2.0, n)
    off  = fill(1.0, n - 1)

    return (1 / h^2) * spdiagm(-1 => off, 0 => main, 1 => off)
end



function dirichlet_laplacian_2d_kron2(Δx, Δy, Nx, Ny)

    # number of interior points
    # Nx = nx - 1
    # Ny = ny - 1

    Dxx = sparse_second_derivative_matrix(Nx, Δx)
    Dyy = sparse_second_derivative_matrix(Ny, Δy)

    Ix = sparse(I, Nx, Nx)
    Iy = sparse(I, Ny, Ny)

    # 2D Laplacian
    L = kron(Iy, Dxx) + kron(Dyy, Ix)

    return L
end

