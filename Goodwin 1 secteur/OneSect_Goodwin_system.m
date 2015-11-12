function [ dX ] = OneSect_Goodwin_system( t,X )

nu = 3;
alpha = 0.025;
beta = 0.02;
delta = 0.01;
phi0 = 0.04/(1-0.04^2);
phi1 = 0.04^3/(1-0.04^2);

K = X(1);
w = X(2);
al = X(3);
N = X(4);

lambda = K/(N*nu/al);
phi_lambda = -phi0 + phi1/(1-lambda)^2;
d_N = beta*N;
d_al = -alpha*al;
d_w = phi_lambda * w;
d_K = (1-w*al)*K/nu - delta*K;

dX = zeros(4,1);
dX(1) = d_K;
dX(2) = d_w;
dX(3) = d_al;
dX(4) = d_N;

end