[T,Z] = euler(@TwoSect_GoodwinKeen_system,[0 100],[0.5 0.5 1.6 1/2 1/5.4*0.675/0.9 0.015 0.015],100000);

nu_1 = 4;
nu_2 = 4;
a_11 = 0.0;
a_12 = 0.0;
a_21 = 0.0;
a_22 = 0.0;

%Q_1 = Z(:,1)/nu_1;
%Q_2 = Z(:,2)/nu_2;
Y_1 = Z(:,1)/nu_1 - a_11*Z(:,1)/nu_1 - a_12*Z(:,2)/nu_2;
Y_2 = Z(:,2)/nu_2 - a_21*Z(:,1)/nu_1 - a_22*Z(:,2)/nu_2;
omega = Z(:,3) .* Z(:,4) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2) ./ (Y_1+Y_2);
lambda = Z(:,4) ./ Z(:,5) .* (Z(:,1)/nu_1 + Z(:,2)/nu_2);
d_1 = Z(:,6)./Y_1;
d_2 = Z(:,7)./Y_2;
d = (Z(:,6)+Z(:,7)) ./ (Y_1+Y_2);

figure
subplot(4,2,1);
plot(T,omega,'-')
legend('\omega')

subplot(4,2,3);
plot(T,lambda,'-')
legend('\lambda')

subplot(4,2,2);
plot(omega,lambda,'-')
legend('\omega \lambda trajectory')

subplot(4,2,4);
plot(T,Z(:,3),'-')
legend('w')

subplot(4,2,5);
plot(T,Z(:,6)./Y_1,'-')
legend('d_1')

subplot(4,2,6);
plot(T,Z(:,7)./Y_2,'-')
legend('d_2')

subplot(4,2,7);
plot(T,d,'-')
legend('d')

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
plot(T,Z(:,6),'-')
legend('D_1')

subplot(5,2,6);
plot(T,Z(:,7),'-')
legend('D_2')

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