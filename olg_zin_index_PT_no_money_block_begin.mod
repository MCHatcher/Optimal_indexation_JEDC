//This mod file was finalised on 9 Feb 2014 (written by M. Hatcher)

var A, b, bd, co, cy, g, k, mu, muy, pi, pro, R, rm, rn, v, q, out, eco, w, rf, rk, btotd, btot, bi, bid, utility, g_ratio, ri, irp, piind;

varexo e, u, u1;
 
parameters alpha, beta, rho, delta, pstar, promean, eps, sharei, N, t, tk;

alpha = 0.263;  //capital share in output
beta = 0.70;    //private discount factor
rho = 0.400;    //persistence in TFP
delta = -15;    //(negative of) coefficient of relative risk aversion
pstar = 0.02;   //yearly inflation target
promean = 1;    //mean of TFP
eps = -0.35;    //intertemporal elasticity of substitution = 1/(1-eps)
//Parameters below are tax rates (t and tk), indexation share (sharei) and a scaling parameter (N) which prevents large numbers
//t and sharei are replaced in each simulation


