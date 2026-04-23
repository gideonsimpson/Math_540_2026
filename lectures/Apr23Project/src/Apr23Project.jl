module Apr23Project

scalar_rhs = (u,p,t) -> u^2 -t;

function lk_rhs!(du, u, p, t)
    du[1]= u[1] -0.1 * u[1]*u[2];
    du[2] = 0.075 * u[1] * u[2] - 1.5 * u[2];
    du
end

export scalar_rhs, lk_rhs!;

end # module Apr23Project
