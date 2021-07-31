//This mod file was finalised on 31 Mar 2014 (written by M. Hatcher)

var A, b, bd, co, cy, g, k, md, m, mu, muy, pi, pro, R, rm, rn, ri, v, q, z, s, x, out, eco, w, rf, rk, btotd, btot, bi, bid, utility, g_ratio, irp, piind;

varexo e, u, u1, u2, u3, u4;
 
parameters alpha, beta, rho, delta, pstar, promean, eps, gamma, sharei, N, t, tk;

alpha = 0.263;  //capital share in output
beta = 0.70;    //private discount factor
rho = 0.400;    //persistence in TFP
delta = -15;    //(negative of) coefficient of relative risk aversion
pstar = 0.02;   //yearly inflation target
promean = 1;    //mean of TFP
eps = -0.35;    //intertemporal elasticity of substitution = 1/(1-eps)
gamma = 0.015;  //real money holdings (cash constraint parameter)
//Parameters below are tax rates (t and tk), indexation share (sharei) and a scaling parameter (N) which prevents large numbers
//t and sharei are replaced in each simulation




