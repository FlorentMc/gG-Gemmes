function [ dX ] = OneSect_GoodwinKeen_system( t,X )

nu_1 = 4;

alpha = 0.025;
beta = 0.02;
delta_1 = 0.01;

theta_1 = 1;

a_11 = 0.0;
a_12 = 0.0;
a_21 = 0.0;
a_22 = 0.0;
phi0 = 0.04/(1-0.04^2);
phi1 = 0.04^3/(1-0.04^2);
r = 0.03;

K_1 = X(1);
K_2 = X(2);
w = X(3);
al = X(4);
N = X(5);
D_1 = X(6);
D_2 = X(7);

Q_1 = K_1/nu_1;

Y_1 = Q_1 - a_11*Q_1;

L_1 = al*Q_1;

Pi_1 = Q_1 - a_11*Q_1  - w*L_1 - r*D_1;


pi = (Pi_1) / (Y_1);
f_pi = -0.0065 + exp(-5)*exp(20*pi);
lambda = al/N*(Q_1);
phi_lambda = -phi0 + phi1/(1-lambda)^2;

d_K_1 = theta_1*f_pi*Y_1 - delta_1*K_1;
d_K_2 = 0;
d_w = phi_lambda * w;
d_al = - alpha * al;
d_N = beta * N;
d_D_1 = theta_1*f_pi*Y_1 - Pi_1;
d_D_2 = 0;

dX = zeros(7,1);
dX(1) = d_K_1;
dX(2) = d_K_2;
dX(3) = d_w;
dX(4) = d_al;
dX(5) = d_N;
dX(6) = d_D_1;
dX(7) = d_D_2;

end