nu_1 = 4;
nu_2 = 4;
alphap = 0.025;
betap = 0.02;
delta_1 = 0.01;
delta_2 = 0.01;
theta_1 = 0.5;
theta_2 = 0.5;
a_11 = 0.01;
a_12 = 0.01;
a_21 = 0.01;
a_22 = 0.01;
phi0 = 0.04/(1-0.04^2);
phi1 = 0.04^3/(1-0.04^2);
r = 0.03;
eta_1 = 4;
eta_2 = 4;
gammainf = 0.8;
zeta_1 = 0.5;
zeta_2 = 0.5;

save parameters.mat nu_1 nu_2 alphap betap delta_1 delta_2 theta_1 theta_2 a_11 a_12 a_21 a_22 phi0 phi1 r eta_1 eta_2 gammainf zeta_1 zeta_2;

periods = 100;
verif_IC=0;
%in order to activate the Initial Conditions verification procedure, uncomment the following line
%verif_IC = 1;
if verif_IC
    numSteps = 2;
else
    numSteps = 1000*periods;
end

[T,Z] = euler(@TwoSect_GoodwinKeenInf_system,[0 periods],[0.3000    0.7000    1.5523    0.5000    0.1389    0.0087    0.0155    1.2000    0.9000],numSteps);

%here, just copy-paste the parameters you chose in the ModelName_system.m
%file

%auxiliary variables
%Q_1 = Z(:,1)/nu_1;
%Q_2 = Z(:,2)/nu_2;
Y_1 = Z(:,1)/nu_1 - a_11*Z(:,1)/nu_1 - a_12*Z(:,2)/nu_2;
Y_2 = Z(:,2)/nu_2 - a_21*Z(:,1)/nu_1 - a_22*Z(:,2)/nu_2;
omega = Z(:,3) .* Z(:,4) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2) ./ (Z(:,8).*Y_1 + Z(:,9).*Y_2);
lambda = Z(:,4) ./ Z(:,5) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2);
d_1 = Z(:,6)./(Z(:,8).*Y_1);
d_2 = Z(:,7)./(Z(:,9).*Y_2);
d = (Z(:,6)+Z(:,7)) ./ (Z(:,8).*Y_1 + Z(:,9).*Y_2);
Pi_1 = Z(:,8).*Z(:,1)/nu_1 - a_11*Z(:,8).*Z(:,1)/nu_1 - a_21*Z(:,9).*Z(:,1)/nu_1 - Z(:,3).*(Z(:,4).*Z(:,1)/nu_1) - r*Z(:,6);
Pi_2 = Z(:,9).*Z(:,2)/nu_2 - a_22*Z(:,9).*Z(:,2)/nu_2 - a_12*Z(:,8).*Z(:,2)/nu_2 - Z(:,3).*(Z(:,4).*Z(:,2)/nu_2) - r*Z(:,7);
pi = (Pi_1+Pi_2) ./ (Z(:,8).*Y_1 + Z(:,9).*Y_2);
inf = zeros(numSteps+1,1);
for k = 1:numSteps
    inf(k) = (zeta_1*(Z(k+1,8)-Z(k,8)) + zeta_2*(Z(k+1,9)-Z(k,9))) ./ (zeta_1*Z(k,8) + zeta_2*Z(k,9));
end
inf(numSteps+1)=inf(numSteps) +inf(numSteps)-inf(numSteps-1);

%procede to verif_IC or plotting procedure
if verif_IC
    omega_init = 0.8; %your wished initial condition for omega
    lambda_init = 0.9; %your wished initial condition for lambda
    d_1_init = 0.1; %your wished initial condition for d_1
    d_2_init = 0.1; %your wished initial condition for d_2
    IC(1) = Z(1,1); %your wished initial condition for K_1
    IC(2) = Z(1,2); %your wished initial condition for K_2
    Q_1_init = IC(1)/nu_1;
    Q_2_init = IC(2)/nu_2;
    IC(8) = Z(1,8); %your wished initial condition for p_1
    IC(9) = Z(1,9); %your wished initial condition for p_2
    IC(4) = 1/2; %arbitrarily set al_init = 1/2
    IC(3) = omega_init*(IC(8)*(1-a_11-a_21)*Q_1_init + IC(9)*(1-a_12-a_22)*Q_2_init) / (IC(4)*(Q_1_init+Q_2_init)); %deduce w_init from other IC
    IC(5) = IC(4)*(Q_1_init+Q_2_init)/lambda_init; %deduce N_init from other IC
    IC(6) = d_1_init*IC(8)*Y_1(1); %deduce D_1_init from other IC
    IC(7) = d_2_init*IC(9)*Y_2(1); %deduce D_2_init from other IC

    disp('IC vector K_1 K_2 w al N D_1 D_2 p_1 p_2'); IC
    
else
    figure
    subplot(6,2,1);
    plot(T,omega,'-')
    legend('\omega')

    subplot(6,2,3);
    plot(T,lambda,'-')
    legend('\lambda')

    subplot(6,2,2);
    plot(omega,lambda,'-')
    legend('\omega \lambda trajectory')

    subplot(6,2,4);
    plot(T,Z(:,3),'-')
    legend('w')

    subplot(6,2,5);
    plot(T,Z(:,6)./Y_1,'-')
    legend('d_1')

    subplot(6,2,6);
    plot(T,Z(:,7)./Y_2,'-')
    legend('d_2')

    subplot(6,2,7);
    plot(T,d,'-')
    legend('d')
    
    subplot(6,2,8);
    plot(T,pi,'-')
    legend('pi')

    subplot(6,2,9);
    plot(T,Pi_1./Y_1,'-')
    legend('pi_1')
    
    subplot(6,2,10);
    plot(T,Pi_2./Y_2,'-')
    legend('pi_2')
    
    subplot(6,2,11);
    plot(T,inf,'-')
    legend('inflation')
    
    subplot(6,2,12);
    plot(T,Z(:,8)./Z(:,9),'-')
    legend('p_1/p_2')
    
    figure
    subplot(7,2,1);
    plot(T,Y_1,'-')
    legend('Y_1')

    subplot(7,2,2);
    plot(T,Y_2,'-')
    legend('Y_2')

    subplot(7,2,3);
    plot(T,Z(:,1),'-')
    legend('K_1')

    subplot(7,2,4);
    plot(T,Z(:,2),'-')
    legend('K_2')

    subplot(7,2,5);
    plot(T,Z(:,6),'-')
    legend('D_1')

    subplot(7,2,6);
    plot(T,Z(:,7),'-')
    legend('D_2')

    subplot(7,2,7);
    plot(T,Z(:,3).*Z(:,4).*(Z(:,1)/nu_1)./Y_1,'-')
    legend('w*L_1/Y_1')

    subplot(7,2,8);
    plot(T,Z(:,3).*Z(:,4).*(Z(:,2)/nu_2)./Y_2,'-')
    legend('w*L_2/Y_2')

    subplot(7,2,9);
    plot(T,Z(:,4).*(Z(:,1)/nu_1)./Z(:,5),'-')
    legend('L_1/N')

    subplot(7,2,10);
    plot(T,Z(:,4).*(Z(:,2)/nu_2)./Z(:,5),'-')
    legend('L_2/N')
    
    subplot(7,2,11);
    plot(T,Pi_1,'-')
    legend('Pi_1')
    
    subplot(7,2,12);
    plot(T,Pi_2,'-')
    legend('Pi_2')
    
    subplot(7,2,13);
    plot(T,Z(:,8),'-')
    legend('p_1')
    
    subplot(7,2,14);
    plot(T,Z(:,9),'-')
    legend('p_2')
    
    figure
    plot3(omega,lambda,d)
    xlabel('\omega')
    ylabel('\lambda')
    zlabel('d')
end