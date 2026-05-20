using SparseArrays

function assemble_laplacian2d(Δx, Δy, nx, ny)
    # Construct the 2D Laplacian for the rectangle with Dirichlet boundary conditions

    rows = Int[] # rows
    cols = Int[] # columuns
    vals = Float64[]  # values

    for k in 1:nx*ny
        # get local coordinates from global coordinate
        j, i = fldmod1(k, nx)
        # corners, then edges, then interior
        if (i == 1 && j == 1)
            append!(rows, k * ones(Int, 3))
            append!(cols, [1, 2, 1 + nx])
            append!(vals, [2 / Δx^2 + 2 / Δy^2, -1 / Δx^2, -1 / Δy^2])
        elseif (i == 1 && j == ny)
            append!(rows, k * ones(Int, 3))
            append!(cols, [1 + (ny - 2) * nx, 1 + (ny - 1) * nx, 2 + (ny - 1) * nx])
            append!(vals, [-1 / Δy^2, 2 / Δx^2 + 2 / Δy^2, -1 / Δx^2])
        elseif (i == nx && j == 1)
            append!(rows, k * ones(Int, 3))
            append!(cols, [nx - 1, nx, 1 + nx])
            append!(vals, [-1 / Δx^2, 2 / Δx^2 + 2 / Δy^2, -1 / Δy^2])
        elseif (i == nx && j == ny)
            append!(rows, k * ones(Int, 3))
            append!(cols, [(nx - 1) * ny, nx * ny - 1, nx * ny])
            append!(vals, [-1 / Δy^2, -1 / Δx^2, 2 / Δx^2 + 2 / Δy^2])
        elseif (i == 1 && j > 1 && j < ny)
            # right edge of domain, (1<j<ny)
            append!(rows, k * ones(Int, 4))
            append!(cols, [1 + (j - 2) * nx, 1 + (j - 1) * nx, 2 + (j - 1) * nx, 1 + j * nx])
            append!(vals, [-1 / Δy^2, 2 / Δx^2 + 2 / Δy^2, -1 / Δx^2, -1 / Δy^2])
        elseif (i == nx && j > 1 && j < ny)
            # left edge of domain, (1<j<ny)
            append!(rows, k * ones(4))
            append!(cols, [(j - 1) * nx, j * nx - 1, j * nx, (j + 1) * nx])
            append!(vals, [-1 / Δy^2, -1 / Δx^2, 2 / Δx^2 + 2 / Δy^2, -1 / Δy^2])
        elseif (i > 1 && i < nx && j == 1)
            # bottom edge of domain, (1<i<nx)
            append!(rows, k * ones(Int, 4))
            append!(cols, [i - 1, i, i + 1, i + nx])
            append!(vals, [-1 / Δx^2, 2 / Δx^2 + 2 / Δy^2, -1 / Δx^2, -1 / Δy^2])
        elseif (i > 1 && i < nx && j == ny)
            # top edge of domain, (1<i<nx)
            append!(rows, k * ones(Int, 4))
            append!(cols, [i + (ny - 2) * nx, i - 1 + (ny - 1) * nx, i + (ny - 1) * nx, i + 1 + (ny - 1) * nx])
            append!(vals, [-1 / Δy^2, -1 / Δx^2, 2 / Δx^2 + 2 / Δy^2, -1 / Δx^2])
        else
            # interior of domain
            append!(rows, (i + (j - 1) * nx) * ones(Int, 5))
            append!(cols, [i + (j - 2) * nx, i - 1 + (j - 1) * nx, i + (j - 1) * nx, i + 1 + (j - 1) * nx, i + j * nx])
            append!(vals, [-1 / Δy^2, -1 / Δx^2, 2 / Δx^2 + 2 / Δy^2, -1 / Δx^2, -1 / Δy^2])
        end
    end

    A = sparse(rows, cols, vals, nx * ny, nx * ny)
    #     display(cols)
    return A
end
