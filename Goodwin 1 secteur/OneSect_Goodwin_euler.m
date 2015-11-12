verif_IC=0;
%in order to activate the Initial Conditions verification procedure, uncomment the following line
%verif_IC = 1;
if verif_IC
    numSteps = 1;
else
    numSteps = 100000;
end

[T,Z] = euler(@OneSect_Goodwin_system,[0 50],[1.0000    1.6000    0.5000    0.1852],numSteps);

%here, just copy-paste the parameters you chose in the ModelName_system.m
%file
nu = 3;
alpha = 0.025;
beta = 0.02;
delta = 0.01;
phi0 = 0.04/(1-0.04^2);
phi1 = 0.04^3/(1-0.04^2);

%auxiliary variables
omega = Z(:,2).*Z(:,3);
lambda = Z(:,1).*Z(:,3)./(nu*Z(:,4));

%procede to verif_IC or plotting procedure
if verif_IC
    omega_init = 0.8; %your wished initial condition for omega
    lambda_init = 0.9; %your wished initial condition for lambda
    IC(1) = 1; %arbitrarily set K_init = 1
    IC(3) = 1/2; %arbitrarily set al_init = 1/2
    IC(2) = omega_init/IC(3); %deduce w_init from other IC
    IC(4) = IC(1)*IC(3)/(nu*lambda_init); %deduce N_init from other IC
    disp('IC vector K w al N'); IC
    
else
    figure
    subplot(2,2,1);
    plot(T,omega,'-')
    legend('\omega')

    subplot(2,2,2);
    plot(T,lambda,'-')
    legend('\lambda')

    subplot(2,2,3);
    plot(omega,lambda,'-')
    legend('\omega \lambda trajectory')
end