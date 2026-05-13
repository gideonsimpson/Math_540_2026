module May12Project

using SparseArrays

function sparse_second_derivative_matrix(n, h)
    @assert n >= 1 "n must be at least 1"
    @assert h > 0 "h must be positive"

    main = fill(-2.0, n)
    off  = fill(1.0, n - 1)

    return (1 / h^2) * spdiagm(-1 => off, 0 => main, 1 => off)
end

export sparse_second_derivative_matrix

end # module May12Project
