close all
rho0=10.963;

figure

T1=linspace(273,923);
rho1 = rho0*(9.9734*10^(-1) - 9.802*10^(-6).*(T1) - 2.705*10^(-10).*(T1.^2) + 4.391*10^(-13).*(T1.^3)).^(-3);
plot(T1,rho1,'b','linewidth',2)
hold on
T2=linspace(923,3120);
rho2= rho0*(9.9672*10^(-1) - 1.179*10^(-5).*(T2) - 2.429*10^(-9).*(T2.^2) + 1.219*10^(-12).*(T2.^3) ).^(-3);
plot(T2,rho2,'r','linewidth',2)
T1=linspace(923,3120);
rho1 = rho0*(9.9734*10^(-1) - 9.802*10^(-6).*(T1) - 2.705*10^(-10).*(T1.^2) + 4.391*10^(-13).*(T1.^3)).^(-3);
plot(T1,rho1,'b--','linewidth',2)

T2=linspace(273,923);
rho2= rho0*(9.9672*10^(-1) - 1.179*10^(-5).*(T2) - 2.429*10^(-9).*(T2.^2) + 1.219*10^(-12).*(T2.^3) ).^(-3);
% plot(T2,rho2,'-.','Color',[0.5 0.5 0.5],'linewidth',2)
plot(T2,rho2,'r--','linewidth',2)




title('Densidad UO2 en función de la temperatura')
xlabel('Temperatura (K)')
ylabel('Densidad UO2 (g/cm3)')
legend('\rho para 273<T<923','\rho para 923<T<3120')
