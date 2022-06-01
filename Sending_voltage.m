syms u(t) v(t);

ode1 = diff(u,2)* 0.1  + 7 * diff(u) + u / (0.0000001) + 0.09 * diff(v,2)   == 7 * cos(2 * pi * 3950*t);
ode2 = diff(v,2)* 0.1  + 10000 * diff(v) + v / (0.0000001) + 0.09 * diff(u,2)   == 0;
odes = [ode1; ode2];
cond1 = u(0) == 0;
cond2 = v(0) == 0;
cond3 = diff(u(0)) == 0;
cond4 = diff(v(0)) == 0;
conds = [cond1; cond2;cond3;cond4];
[uSol(t),vSol(t)] = dsolve(odes)
[uSol(t),vSol(t)] = dsolve(odes,conds)

fplot(uSol)
hold on
fplot(vSol)
grid on
legend('uSol','vSol','Location','best')