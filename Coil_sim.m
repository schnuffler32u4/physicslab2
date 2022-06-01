clc
clearvars
syms a;
syms b;
syms x;
R = 0.105;
u0 = 1.256637061 * 10^(-6);
syms m(a,b,x)
m(a,b,x) == R^2 * cos(a-b) / (sqrt(2*R^2 * (1-cos(a-b)) + x^2));
%M = integral2(m(a,b,x),0,2*pi,0,2*pi)

syms q(t);
ode = diff(q,4) + diff(q,3) + diff(q,2) + diff(q) + q == 5 * cos(55000 * t);
cond1 = q(0) == 0;
cond2 = diff(q(0)) == 0;
cond3 = diff(q(0),2) == 0;
cond4 = diff(q(0),3) == 0;
usol(t) = dsolve(ode)
conds = [cond1; cond2;cond3;cond4];
usoln(t) = dsolve(ode,conds);
fplot(usoln)
grid on