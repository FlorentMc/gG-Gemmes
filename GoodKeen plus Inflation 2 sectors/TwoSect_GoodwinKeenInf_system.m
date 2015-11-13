function [ dX ] = TwoSect_GoodwinKeenInf_system( t,X )

load parameters.mat;

K_1 = X(1);
K_2 = X(2);
w = X(3);
al = X(4);
N = X(5);
D_1 = X(6);
D_2 = X(7);
p_1 = X(8);
p_2 = X(9);

Q_1 = K_1/nu_1;
Q_2 = K_2/nu_2;
Y_1 = Q_1 - a_11*Q_1 - a_12*Q_2;
Y_2 = Q_2 - a_21*Q_1 - a_22*Q_2;
L_1 = al*Q_1;
L_2 = al*Q_2;
Pi_1 = p_1*Q_1 - p_1*a_11*Q_1 - p_2*a_21*Q_1 - w*L_1 - r*D_1;
Pi_2 = p_2*Q_2 - p_2*a_22*Q_2 - p_1*a_12*Q_2 - w*L_2 - r*D_2;
r_1 = (Pi_1 + r*D_1)/(p_1*K_1);
r_2 = (Pi_2 + r*D_2)/(p_1*K_2);
rbar = mean([r_1 r_2]);
CPI = zeta_1*p_1 + zeta_2*p_2;

pi = (Pi_1+Pi_2) / (p_1*Y_1 + p_2*Y_2);
f_pi = -0.0065 + exp(-5)*exp(20*pi);
lambda = al/N*(Q_1+Q_2);
phi_lambda = -phi0 + phi1/(1-lambda)^2;

d_K_1 = theta_1*f_pi*Y_1 - delta_1*K_1;
d_K_2 = theta_2*f_pi*Y_1 - delta_2*K_2;
d_p_1 = eta_1*p_1*nu_1*(rbar-r_1);
d_p_2 = eta_2*p_1*nu_2*(rbar-r_2);

inf = (zeta_1*d_p_1 + zeta_2*d_p_2)/CPI;

d_w = (phi_lambda + gammainf*inf) * w;
d_al = - alphap * al;
d_N = betap * N;
d_D_1 = theta_1*f_pi*Y_1 - Pi_1;
d_D_2 = theta_2*f_pi*Y_1 - Pi_2;

dX = zeros(9,1);
dX(1) = d_K_1;
dX(2) = d_K_2;
dX(3) = d_w;
dX(4) = d_al;
dX(5) = d_N;
dX(6) = d_D_1;
dX(7) = d_D_2;
dX(8) = d_p_1;
dX(9) = d_p_2;

end