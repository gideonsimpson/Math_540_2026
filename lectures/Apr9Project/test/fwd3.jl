f = (t,u) -> u;

u0 = 1.0
t0 = 0.0;
tmax = 1.0;
N = 10;
t_vals, u_vals = explicit_euler(f, u0, t0, tmax, N);
u_exact = (1+1/N)^N;
isapprox(u_vals[end], u_exact)