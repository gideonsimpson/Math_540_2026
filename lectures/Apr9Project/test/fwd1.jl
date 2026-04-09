# f(t,u) = 0.0
f = (t,u) -> 0.0;
u0 = 1.0
t0 = 0.0;
tmax = 1.0;
N = 10;
t_vals, u_vals = explicit_euler(f, u0, t0, tmax, N);
# check that u(tmax)=1 is true:
isapprox(u_vals[end], 1.0)