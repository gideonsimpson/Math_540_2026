f = (t,u) -> u;

u0 = 1.0
t0 = 0.0;
tmax = 1.0;
N = 10^6;
t_vals, u_vals = explicit_euler(f, u0, t0, tmax, N);
# check error is smaller than 10^-4, for h = 10^-6
abs(u_vals[end] - exp(1)) < 1e-4
