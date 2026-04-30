nx = 100;
x = LinRange(-5, 5, nx + 1)[2:end]; # skip the first point because we are not solving there
m = 3.0; # slope
u = m * (x .- (-5));
dx = x[2] - x[1];
Dx = sparse_backwards_difference_matrix(nx, dx);

# compare against the exact discrete derivative, m * 1 vector
norm(Dx * u - m * ones(nx),Inf) <  1e-12
