%XSTEAM
format long
T=[];
close all
clear vars
c=0;
ht=zeros(1,length(T));
rho=zeros(1,length(T));
mu=zeros(1,length(T));
lambda=zeros(1,length(T));
cp=zeros(1,length(T));
G_agua_total  = 12893; %kg/s
n_fa=221;

p=12.655e-3;
R_agua = sqrt(((p^2))/pi);
R_clad_ext = 9.517*0.001/2;
S_agua= 264*221*pi*(R_agua^2 - R_clad_ext^2);

u=G_agua_total/(S_agua*XSteam('rho_pt',155,286));

Str = p^2 - pi*R_clad_ext^2;
Ptr = 4*p - 2*pi*R_clad_ext;
Dh=4*Str/Ptr;

for Ti=286:0.5:floor(XSteam('Tsat_p',155))
    T=[T,Ti];
    c=c+1;
    ht(c)=XSteam('h_pt',155,Ti);
    rho(c)=XSteam('rho_pt',155,Ti);
    mu(c)=XSteam('my_pt',155,Ti);
    lambda(c) = XSteam('tc_pt',155,Ti);
    cp(c)=XSteam('Cp_pt',155,Ti)*1000;
end
nu = mu./rho;

Re = u.*Dh./nu;
Pr = mu.*cp./lambda;

Nu=0.023.*(Re.^(0.8)).*(Pr.^(0.4));
h = lambda.*Nu / Dh;

Tk=T+273.15;
h_p=polyfit(Tk,ht,2);
h_val=polyval(h_p,Tk);

rho_p=polyfit(Tk,rho,2);
rho_val=polyval(rho_p,Tk);

%PLOTS
% figure
% plot(Tk(1:10:end),ht(1:10:end),'r+','markersize',8,'linewidth',1.2)
% hold on
% plot(Tk,h_val,'b','linewidth',1.5)
% title('Entalpía del agua en función de la temperatura a 155bar')
% legend('XSteam','Polynomial')
% xlabel('Temperatura (K)')
% ylabel('Entalpía (kJ/kg)')
% legend('Datos XSteam','Polinomio')
% grid on
% display(100*max(abs(h-h_val)./h))
% 
% %PLOTS
% figure
% grid on
% plot(Tk(1:10:end),rho(1:10:end),'r+','markersize',8,'linewidth',1.2)
% hold on
% plot(Tk,rho_val,'b','linewidth',1.5)
% title('Densidad del agua en función de la temperatura a 155bar')
% legend('XSteam','Polynomial')
% xlabel('Temperatura (K)')
% ylabel('Densidad (kg/m^3)')
% display(100*max(abs(rho-rho_val)./rho))
% legend('Datos XSteam','Polinomio')
% grid on
% figure
% grid on
% plot(Tk,Nu,'bx','linewidth',1.2)
% xlabel('Temperatura (K)')
% ylabel('Nu')
% title('Número de Nusselt en función de la temperatura')
% grid on
% figure
% grid on
% plot(Tk,Pr,'bx','linewidth',1.2)
% xlabel('Temperatura (K)')
% ylabel('Pr')
% title('Número de Prandtl en función de la temperatura')
% grid on
% figure
% grid on
% plot(Tk,Re,'bx','linewidth',1.2)
% xlabel('Temperatura (K)')
% ylabel('Re')
% title('Número de Reynolds en función de la temperatura')
% grid on
% figure
% 
% plot(Tk,h,'bx','linewidth',1.2)
% xlabel('Temperatura (K)')
% ylabel('h (W/Km2)')
% title('Coeficiente de convección en función de la temperatura')
% grid on
