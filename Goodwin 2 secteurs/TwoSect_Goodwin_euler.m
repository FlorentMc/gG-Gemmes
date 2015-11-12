verif_IC=0;
%in order to activate the Initial Conditions verification procedure, uncomment the following line
%verif_IC = 1;
if verif_IC
    numSteps = 1;
else
    numSteps = 100000;
end

[T,Z] = euler(@TwoSect_Goodwin_system,[0 100],[0.3000    0.7000    1.5680    0.5000    0.1713],numSteps);

%here, just copy-paste the parameters you chose in the ModelName_system.m
%file
nu_1 = 3;
nu_2 = 3;
alpha = 0.025;
beta = 0.02;
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

%auxiliary variables
%Q_1 = Z(:,1)/nu_1;
%Q_2 = Z(:,2)/nu_2;
Y_1 = Z(:,1)/nu_1 - a_11*Z(:,1)/nu_1 - a_12*Z(:,2)/nu_2;
Y_2 = Z(:,2)/nu_2 - a_21*Z(:,1)/nu_1 - a_22*Z(:,2)/nu_2;
omega = Z(:,3) .* Z(:,4) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2) ./ (Y_1+Y_2);
lambda = Z(:,4) ./ Z(:,5) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2);

%procede to verif_IC or plotting procedure
if verif_IC
    omega_init = 0.8; %your wished initial condition for omega
    lambda_init = 0.9; %your wished initial condition for lambda
    IC(1) = Z(1,1); %your wished initial condition for K_1
    IC(2) = Z(1,2); %your wished initial condition for K_2
    Q_1_init = IC(1)/nu_1;
    Q_2_init = IC(2)/nu_2;
    IC(4) = 1/2; %arbitrarily set al_init = 1/2
    IC(3) = omega_init*((1-a_11-a_21)*Q_1_init + (1-a_12-a_22)*Q_2_init) / (IC(4)*(Q_1_init+Q_2_init)); %deduce w_init from other IC
    IC(5) = IC(4)*(Q_1_init+Q_2_init)/lambda_init; %deduce N_init from other IC
    disp('IC vector K_1 K_2 w al N'); IC
    
else
    figure
    subplot(2,2,1);
    plot(T,omega,'-')
    legend('\omega')

    subplot(2,2,3);
    plot(T,lambda,'-')
    legend('\lambda')

    subplot(2,2,2);
    plot(omega,lambda,'-')
    legend('\omega \lambda trajectory')

    figure
    subplot(5,2,1);
    plot(T,Y_1,'-')
    legend('Y_1')

    subplot(5,2,2);
    plot(T,Y_2,'-')
    legend('Y_2')

    subplot(5,2,3);
    plot(T,Z(:,1),'-')
    legend('K_1')

    subplot(5,2,4);
    plot(T,Z(:,2),'-')
    legend('K_2')

    subplot(5,2,5);
    plot(T,Z(:,3),'-')
    legend('w')

    subplot(5,2,7);
    plot(T,Z(:,3).*Z(:,4).*(Z(:,1)/nu_1)./Y_1,'-')
    legend('w*L_1/Y_1')

    subplot(5,2,8);
    plot(T,Z(:,3).*Z(:,4).*(Z(:,2)/nu_2)./Y_2,'-')
    legend('w*L_2/Y_2')

    subplot(5,2,9);
    plot(T,Z(:,4).*(Z(:,1)/nu_1)./Z(:,5),'-')
    legend('L_1/N')

    subplot(5,2,10);
    plot(T,Z(:,4).*(Z(:,2)/nu_2)./Z(:,5),'-')
    legend('L_2/N')
end