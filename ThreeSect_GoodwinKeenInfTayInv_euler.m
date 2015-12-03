periods = 100;
verif_IC=0;
%in order to activate the Initial Conditions verification procedure, uncomment the following line
%verif_IC = 1;
if verif_IC
    numSteps = 2;
else
    numSteps = 1000*periods;
end

[T,Z] = euler(@ThreeSect_GoodwinKeenInfTayInv_system,[0 periods],[1/3  1/3  1/3  1.6000  0.5000  0.1389  0.0083  0.0083  0.0083   1  1  1  1/3  1/3  1/3],numSteps);

%here, just copy-paste the parameters you chose in the ModelName_system.m
%file
nu_1 = 4;
nu_2 = 4;
nu_3 = 4;
alpha = 0.025;
beta = 0.02;
delta_1 = 0.01;
delta_2 = 0.01;
delta_3 = 0.01;
a_11 = 0.0;
a_12 = 0.0;
a_13 = 0.0;
a_21 = 0.0;
a_22 = 0.0;
a_23 = 0.0;
a_31 = 0.0;
a_32 = 0.0;
a_33 = 0.0;
phi0 = 0.04/(1-0.04^2);
phi1 = 0.04^3/(1-0.04^2);
eta_1 = 0.1;
eta_2 = 0.1;
eta_3 = 0.1;
gamma = 0.8;
zeta_1 = 1/3;
zeta_2 = 1/3;
zeta_3 = 1/3;
rstar = 0.03; %"neutral" short term interest rate
istar = 0.005; %target inflation rate
u = 9.8813e-324; %auxiliary parameter to make the Taylor rule differentiable
phi_T = 1.5; %Taylor rule reactivity
sigma = 0.1; %adjustment speed parameter for allocation of investment

%auxiliary variables
%Q_1 = Z(:,1)/nu_1;
%Q_2 = Z(:,2)/nu_2;
%Q_3 = Z(:,3)/nu_3;
Y_1 = Z(:,1)/nu_1 - a_11*Z(:,1)/nu_1 - a_12*Z(:,2)/nu_2 - a_13*Z(:,3)/nu_3;
Y_2 = Z(:,2)/nu_2 - a_21*Z(:,1)/nu_1 - a_22*Z(:,2)/nu_2 - a_23*Z(:,3)/nu_3;
Y_3 = Z(:,3)/nu_3 - a_31*Z(:,1)/nu_1 - a_32*Z(:,2)/nu_2 - a_33*Z(:,3)/nu_3;
omega = Z(:,4) .* Z(:,5) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2 + Z(:,3)/nu_3) ./ (Z(:,10).*Y_1+Z(:,11).*Y_2+Z(:,12).*Y_3);
lambda = Z(:,5) ./ Z(:,6) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2 + Z(:,3)/nu_3);
d_1 = Z(:,7)./(Z(:,10).*Y_1);
d_2 = Z(:,8)./(Z(:,11).*Y_2);
d_3 = Z(:,9)./(Z(:,12).*Y_3);
d = (Z(:,7)+Z(:,8)+Z(:,9)) ./ (Z(:,10).*Y_1+Z(:,11).*Y_2+Z(:,12).*Y_3);

inf = zeros(numSteps+1,1);
for k = 1:numSteps
    inf(k) = ((zeta_1*(Z(k+1,10)-Z(k,10)) + zeta_2*(Z(k+1,11)-Z(k,11))+ zeta_3*(Z(k+1,12)-Z(k,12)))/(periods/numSteps)) ./ (zeta_1*Z(k,10) + zeta_2*Z(k,11) + zeta_3*Z(k,12));
end
inf(numSteps+1)=inf(numSteps) +inf(numSteps)-inf(numSteps-1);

iota = rstar + inf + phi_T*(inf-istar);
r = (iota + sqrt(iota.^2 + u))/2;

Pi_1 = Z(:,10).*Z(:,1)/nu_1 - a_11*Z(:,10).*Z(:,1)/nu_1 - a_21*Z(:,11).*Z(:,1)/nu_1  - a_31*Z(:,12).*Z(:,1)/nu_1 - Z(:,4).*(Z(:,5).*Z(:,1)/nu_1) - r.*Z(:,7);
Pi_2 = Z(:,11).*Z(:,2)/nu_2 - a_12*Z(:,10).*Z(:,2)/nu_2 - a_22*Z(:,11).*Z(:,2)/nu_2  - a_32*Z(:,12).*Z(:,2)/nu_2 - Z(:,4).*(Z(:,5).*Z(:,2)/nu_2) - r.*Z(:,8);
Pi_3 = Z(:,12).*Z(:,3)/nu_3 - a_13*Z(:,10).*Z(:,3)/nu_3 - a_23*Z(:,11).*Z(:,3)/nu_3  - a_33*Z(:,12).*Z(:,3)/nu_3 - Z(:,4).*(Z(:,5).*Z(:,3)/nu_3) - r.*Z(:,9);

pi = (Pi_1+Pi_2+Pi_3) ./ (Z(:,10).*Y_1+Z(:,11).*Y_2+Z(:,12).*Y_3);

%procede to verif_IC or plotting procedure
if verif_IC
    omega_init = 0.8; %your wished initial condition for omega
    lambda_init = 0.9; %your wished initial condition for lambda
    d_1_init = 0.1; %your wished initial condition for d_1
    d_2_init = 0.1; %your wished initial condition for d_2
    d_3_init = 0.1; %your wished initial condition for d_3
    IC(1) = Z(1,1); %your wished initial condition for K_1
    IC(2) = Z(1,2); %your wished initial condition for K_2
    IC(3) = Z(1,3); %your wished initial condition for K_3
    Q_1_init = IC(1)/nu_1;
    Q_2_init = IC(2)/nu_2;
    Q_3_init = IC(3)/nu_3;
    IC(10) = Z(1,10); %your wished initial condition for p_1
    IC(11) = Z(1,11); %your wished initial condition for p_2
    IC(12) = Z(1,12); %your wished initial condition for p_3
    IC(13) = Z(1,13); %your wished initial condition for theta_1
    IC(14) = Z(1,14); %your wished initial condition for theta_2
    IC(15) = Z(1,15); %your wished initial condition for theta_2
    IC(5) = 1/2; %arbitrarily set al_init = 1/2
    IC(4) = omega_init*(IC(10)*Y_1(1)+IC(11)*Y_2(1)+IC(12)*Y_3(1)) / (IC(5)*(Q_1_init+Q_2_init+Q_3_init)); %deduce w_init from other IC
    IC(6) = IC(5)*(Q_1_init+Q_2_init+Q_3_init)/lambda_init; %deduce N_init from other IC
    IC(7) = d_1_init*IC(10)*Y_1(1); %deduce D_1_init from other IC
    IC(8) = d_2_init*IC(11)*Y_2(1); %deduce D_2_init from other IC
    IC(9) = d_3_init*IC(12)*Y_3(1); %deduce D_3_init from other IC

    disp('IC vector K_1 K_2 K_3 w al N D_1 D_2 D_3 p_1 p_2 p_3 theta_1 theta_2 theta_3'); IC
    
else
    figure
    subplot(4,3,1);
    plot(T,omega,'-')
    legend('\omega')
    
    subplot(4,3,2);
    plot(T,lambda,'-')
    legend('\lambda')
        
    subplot(4,3,3);
    plot(omega,lambda,'-')
    legend('\omega \lambda trajectory')

    subplot(4,3,4);
    plot(T,Z(:,4),'-')
    legend('w')
        
    subplot(4,3,5);
    plot(T,d,'-')
    legend('d')
    
    subplot(4,3,6);
    plot(T,pi,'-')
    legend('\pi')
    
    subplot(4,3,7);
    plot(T,d_1,'-')
    legend('d_1')

    subplot(4,3,8);
    plot(T,d_2,'-')
    legend('d_2')
    
    subplot(4,3,9);
    plot(T,d_3,'-')
    legend('d_3')

    subplot(4,3,10);
    plot(T,Pi_1./(Z(:,10).*Y_1),'-')
    legend('\pi_1')
    
    subplot(4,3,11);
    plot(T,Pi_2./(Z(:,11).*Y_2),'-')
    legend('\pi_2')
    
    subplot(4,3,12);
    plot(T,Pi_3./(Z(:,12).*Y_3),'-')
    legend('\pi_3')
   
    
    figure
    subplot(4,3,1);
    plot(T,Y_1,'-')
    legend('Y_1')

    subplot(4,3,2);
    plot(T,Y_2,'-')
    legend('Y_2')
    
    subplot(4,3,3);
    plot(T,Y_3,'-')
    legend('Y_3')

    subplot(4,3,4);
    plot(T,Z(:,1),'-')
    legend('K_1')

    subplot(4,3,5);
    plot(T,Z(:,2),'-')
    legend('K_2')
    
    subplot(4,3,6);
    plot(T,Z(:,3),'-')
    legend('K_3')

    subplot(4,3,7);
    plot(T,Z(:,7),'-')
    legend('D_1')

    subplot(4,3,8);
    plot(T,Z(:,8),'-')
    legend('D_2')
    
    subplot(4,3,9);
    plot(T,Z(:,9),'-')
    legend('D_3')
    
    subplot(4,3,10);
    plot(T,Pi_1,'-')
    legend('Pi_1')
    
    subplot(4,3,11);
    plot(T,Pi_2,'-')
    legend('Pi_2')
    
    subplot(4,3,12);
    plot(T,Pi_3,'-')
    legend('Pi_3')

    
    figure
    subplot(5,3,1);
    plot(T,Z(:,10),'-')
    legend('p_1')
    
    subplot(5,3,2);
    plot(T,Z(:,11),'-')
    legend('p_2')
    
    subplot(5,3,3);
    plot(T,Z(:,12),'-')
    legend('p_3')
    
    subplot(5,3,4);
    plot(T,inf,'-')
    legend('inflation')
    
    subplot(5,3,5);
    plot(T,Z(:,11)./Z(:,10),'-')
    legend('p_2/p_1')
    
    subplot(5,3,6);
    plot(T,Z(:,12)./Z(:,10),'-')
    legend('p_3/p_1')
    
    subplot(5,3,7);
    plot(T,r,'-')
    legend('r')
    
    subplot(5,3,8);
    plot(T,iota,'-')
    legend('iota')
    
    subplot(5,3,10);
    plot(T,Z(:,4).*Z(:,5).*(Z(:,1)/nu_1)./(Z(:,10).*Y_1),'-')
    legend('\omega_1')
    
    subplot(5,3,11);
    plot(T,Z(:,4).*Z(:,5).*(Z(:,2)/nu_2)./(Z(:,11).*Y_2),'-')
    legend('\omega_2')
    
    subplot(5,3,12);
    plot(T,Z(:,4).*Z(:,5).*(Z(:,3)/nu_3)./(Z(:,12).*Y_3),'-')
    legend('\omega_3')

    subplot(5,3,13);
    plot(T,(Z(:,1)/nu_1)./(Z(:,1)/nu_1 + Z(:,2)/nu_2 + Z(:,3)/nu_3),'-')
    legend('L_1/L')

    subplot(5,3,14);
    plot(T,(Z(:,2)/nu_2)./(Z(:,1)/nu_1 + Z(:,2)/nu_2 + Z(:,3)/nu_3),'-')
    legend('L_2/L')
    
    subplot(5,3,15);
    plot(T,(Z(:,3)/nu_3)./(Z(:,1)/nu_1 + Z(:,2)/nu_2 + Z(:,3)/nu_3),'-')
    legend('L_3/L')
    
    
    figure
    subplot(2,3,1);
    plot(T,(Z(:,10).*Z(:,1)/nu_1 - a_11*Z(:,10).*Z(:,1)/nu_1 - a_21*Z(:,11).*Z(:,1)/nu_1 - a_31*Z(:,12).*Z(:,1)/nu_1 - Z(:,4).*(Z(:,5).*Z(:,1)/nu_1))./(Z(:,10).*Z(:,1)),'-')
    legend('r_1')
    
    subplot(2,3,2);
    plot(T,(Z(:,11).*Z(:,2)/nu_2 - a_12*Z(:,10).*Z(:,2)/nu_2 - a_22*Z(:,11).*Z(:,2)/nu_2 - a_32*Z(:,12).*Z(:,2)/nu_2 - Z(:,4).*(Z(:,5).*Z(:,2)/nu_2))./(Z(:,10).*Z(:,2)),'-')
    legend('r_2')
    
    subplot(2,3,3);
    plot(T,(Z(:,12).*Z(:,3)/nu_3 - a_13*Z(:,10).*Z(:,3)/nu_3 - a_23*Z(:,11).*Z(:,3)/nu_3 - a_33*Z(:,12).*Z(:,3)/nu_3 - Z(:,4).*(Z(:,5).*Z(:,3)/nu_3))./(Z(:,10).*Z(:,3)),'-')
    legend('r_3')
    
    subplot(2,3,4);
    plot(T,Z(:,13),'-')
    legend('\theta_1')
    
    subplot(2,3,5);
    plot(T,Z(:,14),'-')
    legend('\theta_2')
    
    subplot(2,3,6);
    plot(T,Z(:,15),'-')
    legend('\theta_3')
    
    figure
    plot3(omega,lambda,d)
    xlabel('\omega')
    ylabel('\lambda')
    zlabel('d')
end
