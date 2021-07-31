tk = 2.3*t; 
N = 100000;

model;

//Budget constraint of young consumers
cy = w*(1-t) - btotd - k - md; 

//Budget constraint of old consumers
co = (1-tk)*A*k(-1) + ( sharei*rf(-1)*rm*piind + (1-sharei)*R(-1)*rm )*btotd(-1) + rm*md(-1);

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
muy = mu(+1);

//Money demand/CIA constraint
md = gamma;

//Money market equilibrium
m = md;

//Inflation shock
v = u;

//Shock to indexed inflation
q = u1;

//Money supply rule
//Horizon of 1 year
m = m(-1)*(1+pi)^(-1) * ( (1+pstar)^20 )*(1+v)/(1+v(-1));

//Indexed inflation
piind =( (1+pstar)^20 )*(1+q)/(1+q(-1));

//Productivity (ie TFP)
pro = exp(e)*(pro(-1)^rho)*(promean^(1-rho));

//Euler equation for capital
muy = beta*(1-tk)*mu(+1)*pro(+1)*alpha*k^(alpha-1);

//Real return on money balances
rm = 1/(1+pi);

//Euler equation for nominal bonds 
muy = beta*R*mu(+1)*rm(+1);

//Euler equation for indexed bonds
muy = beta*rf*mu(+1)*rm(+1)*piind(+1);

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
mu = (1/N)*( co^delta )/ ( (eco(-1) / (10^(1+delta)) )^((1-eps+delta)/(1+delta)) );

//Expected consumption term in utility
//Scaling by 10 prevents large numbers
eco = (10*co(+1))^(1+delta);

//Epstein-Zin preferences
//Scaling by 1/(10^30) mean that initial value for utility does not need to changed with small changes on calibration 
//The same scaling is applied under IT
utility = (1/(10^30))*(1/(1+delta))*( (cy)^(eps) + beta*(eco/(10^(1+delta)))^(eps/(1+delta)) )^((1+delta)/eps);

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
cy = 0.1842; 
R = 2.5877;
co = cy;
muy = 0.0000982;
mu = muy; 
pi = 0.8114;
pro = promean;
rm = 0.5521;
md = gamma;
m = md;
A = 1.9195;
g = 0.0563;
utility = -0.000000000023;
out = 0.4920;
eco = 0.000194;
w = 0.3626;
rf = rn;
ri = rn;
rk = 0;
g_ratio = 0.1146;
irp = 0;
v = 0;
q = 0;
piind = 1+pi;
t = 0.11;
tk = 2.3*t;
end;



