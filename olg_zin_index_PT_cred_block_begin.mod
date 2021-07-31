//This mod file was finalised on 4 March 2014 (written by M. Hatcher)

var A, b, bd, co, co_IT, cy, g, k, mu, mu_IT, muy, pi, pi_IT, pro, R, rm, rm_IT, rn, v, z, out, eco, eco_IT, w, rf, rk, ri, btotd, btot, bi, bid, utility, g_ratio, irp, piind, piind_IT;

varexo e, x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9, x_910, x_911, x_912, x_913, x_914, x_915, x_916, x_917, x_918, x_919, x_920;

parameters alpha, beta, rho, delta, pstar, promean, eps, m, cred_PLT, sharei, N, t, t_IT, tk, tk_IT;

alpha = 0.263;  //capital share in output
beta = 0.70;    //private discount factor
rho = 0.400;    //persistence in TFP
delta = -15;    //(negative of) coefficient of relative risk aversion
pstar = 0.02;   //yearly inflation target
promean = 1;    //mean of TFP
eps = -0.35;    //intertemporal elasticity of substitution = 1/(1-eps)
m = 0.015;      //real money holdings (cash constraint parameter)
cred_PLT = 0.9; //perceived probability of reversion to IT regime

//Parameters below are tax rates (t and tk), indexation share (sharei) and a scaling parameter (N) which prevents large numbers
//t and sharei are replaced in each simulation





