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





t = 1.116385e-001; 
t_IT = 1.115859e-001; 
sharei = 1; 
tk_IT = 2.3*t_IT;
tk = 2.3*t; 
N = 100000;

model;

//Budget constraint of young consumers
cy = w*(1-t) - btotd - k - m; 

//Budget constraint of old consumers (under PT)
co = (1-tk)*A*k(-1) + ( sharei*rf(-1)*rm*piind + (1-sharei)*R(-1)*rm )*btotd(-1) + rm*m(-1);

//Budget constraint of old under IT
co_IT = (1-tk_IT)*A*k(-1) + ( sharei*rf(-1)*rm_IT*piind_IT + (1-sharei)*R(-1)*rm_IT )*btotd(-1) + rm_IT*m(-1);

//Government budget constraint
g = w*t + tk*A*k(-1) + b - R(-1)*rm*b(-1) + bi - rf(-1)*rm*piind*bi(-1) + m - rm*m(-1);

//Government spending to output ratio
g_ratio = g/out;

//Total demand for bonds
btotd = bd + bid;

//Nominal bonds market clearing
bd = b;

//Indexed bonds market clearing
bid = bi;

//Supply of indexed bonds
bi = sharei*btot;

//Supply of nominal bonds
b = (1-sharei)*btot;

//Total bond supply equation (ie eq for btot) 
//NB This equation assumes that the bond supply is set with knowledge of private agents' perceptions about monetary policy   
muy = (1-cred_PLT)*mu_IT(+1) + cred_PLT*mu(+1);

//Inflation shock
v = x_919;

//Inflation shock 1
z = x_920;

//Money supply rule under PT
m = m(-1)*(1+pi)^(-1) * ( (1+pstar)^20 ) * (1+z)/(1+z(-1));

//Money supply rule under IT
m = m(-1)*(1+pi_IT)^(-1) * ( (1+pstar)^20 )*(1+x_1)*(1+x_2)*(1+x_3)*(1+x_4)*(1+x_5)*(1+x_6)*(1+x_7)*(1+x_8)*(1+x_9)*(1+x_910)*(1+x_911)*(1+x_912)*(1+x_913)*(1+x_914)*(1+x_915)*(1+x_916)*(1+x_917)*(1+x_918)*(1+x_919)*(1+x_920);

//Indexed inflation under IT
piind_IT =( (1+pstar)^20 )*(1+z(-1))*(1+x_1)*(1+x_2)*(1+x_3)*(1+x_4)*(1+x_5)*(1+x_6)*(1+x_7)*(1+x_8)*(1+x_9)*(1+x_910)*(1+x_911)*(1+x_912)*(1+x_913)*(1+x_914)*(1+x_915)*(1+x_916)*(1+x_917)*(1+x_918)*(1+x_919);

//Indexed inflation under PT
piind = ( (1+pstar)^20 ) * (1+v)/(1+v(-1));

//Productivity 
pro = exp(e)*(pro(-1)^rho)*(promean^(1-rho));

//Euler equation for capital
muy = beta*alpha*k^(alpha-1)*((1-cred_PLT)*(1-tk_IT)*pro(+1)*mu_IT(+1) + cred_PLT*(1-tk)*pro(+1)*mu(+1));

//Real return on money balances
rm = 1/(1+pi);

//Real return on money balances under IT
rm_IT = 1/(1+pi_IT);

//Euler equation for nominal bonds 
muy = beta*R*( (1-cred_PLT)*mu_IT(+1)*rm_IT(+1) + cred_PLT*mu(+1)*rm(+1) );

//Euler equation for indexed bonds
muy = beta*rf*( (1-cred_PLT)*mu_IT(+1)*rm_IT(+1)*piind_IT(+1) + cred_PLT*mu(+1)*rm(+1)*piind(+1) );

//Real return on nominal bonds
rn = R(-1)*rm;

//Real return on indexed bonds
ri = rf(-1)*rm*piind;

//Real return on capital
A = pro*alpha*k(-1)^(alpha-1);

//Equation to calculate average annual after-tax risk premium on capital
rk = ( ((1-tk)*A)^(1/20) - 1 ) - (ri^(1/20) - 1);

//Equation to calculate the average diff between the inflation risk premium on nominal and indexed debt
irp = rn - ri;

//LHS term of Euler equation
//Scaling by 1/N prevents large numbers
muy = (1/N)*cy^(eps-1);

//RHS term of Euler equation
//Scaled and unscaled by 10^(1+delta) due to definition of eco (see below)
//Scaling by 1/N prevents large numbers
mu = (1/N)*( co^delta )/ ( ( (1-cred_PLT)*eco_IT(-1) / (10^(1+delta)) + cred_PLT*eco(-1) /(10^(1+delta)) )^((1-eps+delta)/(1+delta)) );

//RHS term of Euler equation under IT
//Scaled and unscaled by 10^(1+delta) due to definition of eco (see below)
//Scaling by 1/N prevents large numbers
mu_IT = (1/N)*( co_IT^delta )/ ( ( (1-cred_PLT)*eco_IT(-1) / (10^(1+delta)) + cred_PLT*eco(-1) /(10^(1+delta)) )^((1-eps+delta)/(1+delta)) );

//Expected consumption term in the sdf
//Scaling by 10 prevents large numbers
eco = (10*co(+1))^(1+delta);

//Expected consumption term in the sdf (under IT)
//Scaling by 10 prevents large numbers
eco_IT = (10*co_IT(+1))^(1+delta);

//Epstein-Zin preferences
//Scaling by 1/(10^30) means that initval for utility does not need to changed with small changes on calibration 
//The same scaling is applied under both regimes
utility = (1/(10^30))*(1/(1+delta))*( (cy)^(eps) + beta*( (1-cred_PLT)*eco_IT/(10^(1+delta)) + cred_PLT*eco/(10^(1+delta)) )^(eps/(1+delta)) )^((1+delta)/eps); 

//Aggregate output
out = pro*k(-1)^alpha; 

//Wages
w = out - A*k(-1);

end;

initval;
btotd = 0.0557;
btot = btotd;
bid = 0.0111;
bd = 0.0446;
bi = bid;
b = bd;
k = 0.0674;
rn = 1/beta;
ri = rn;
cy = 0.1842; 
R = 2.5877;
co = cy;
co_IT = cy;
muy = 0.0000982;
mu = muy;
mu_IT = muy; 
pi = 0.8114;
pi_IT = pi;
pro = promean;
rm = 0.5521;
rm_IT = rm;
A = 1.9195;
g = 0.0563;
utility = -0.000000000023;
out = 0.4920;
eco = 0.000194;
eco_IT = eco;
w = 0.3626;
//lambda = 0.0000602;
rf = rn;
rk = 0;
g_ratio = 0.1146;
irp = 0;
v = 0;
z = 0;
piind = 1+pi;
piind_IT = piind;
end;vcov = [3.100000e-003 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004]; 
