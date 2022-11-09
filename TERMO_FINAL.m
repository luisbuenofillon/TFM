clearvars all
close all
clc
format long
warning('OFF','ALL')
%% Datos geométricos

R_UO2 = 8.239*0.001/2; %m
e=0.571*0.001; %m
R_clad_ext = 9.517*0.001/2; %m
R_clad_int = R_clad_ext - e; %m
n_pins = 264; %uds
R_guide_tubes=12.259*0.001/2; %m
n_guide_tubes = 25; %uds
n_cells =221; %uds
length_cell = 21.606*0.01; %m
H=3.673;

%% Datos termohidráulicos

T_inlet = 286 + 273.15; % K
Pressure = 155*10e5 ;%bar
G_agua_total  = 12893; %kg/s

%% R_agua

S_cell = (length_cell)^2;
S_occ= n_pins*pi()*(R_clad_ext)^2 + n_guide_tubes*pi()*(R_guide_tubes)^2;
S_agua = S_cell - S_occ;
S_pin = pi()*R_clad_ext;
p=12.655e-3;
R_agua = sqrt(((p^2))/pi);

%% Datos iniciales

Tri=900; %K

Potencia=2775e6/ (n_cells*n_pins) ;%W
Vreactor =pi*H*R_clad_ext^2; %m3
pgen = Potencia/Vreactor;
hconv=2.399046893418894e+04;


%% Inicialización

nr=100; % número de puntos perfil radial
dr = R_clad_ext/(nr+1); % paso axial
r = zeros(nr+1,1);
Xr = zeros(nr+1,1);
for k=1:nr+1
    r(k) = dr*k;
    Xr(k) = Tri;
end
Xr(2)=Xr(1);
w0=3.1;

nz=100; % número de puntos perfil axial
dz = H/(nz+1); %paso axial
Xz = zeros(2*nz+2,1);
z = zeros(nz+1,1);
for k=1:nz+1
    z(k) = dz*k;
end
for k=1:2*nz+2
    if k<=nz+1
        Xz(k)=0;
    else
        Xz(k)=w0;
    end    
end

Tam=T_inlet;
globalerr=[];
globalerrz1=[];
globalerrz2=[];

%% Volúmen de control V1
for iteracion=1:20
    %fprintf('Iteración %d/%d....\n',iteracion,20)
    for it=1:10
        Ar=sparse(nr+1,nr+1);
        br = -pgen*ones(nr+1,1);
        br(1)=0;
        br(end)= -hconv*Tam;
        for i=2:nr
            ri=r(i);
            Tp10 = Xr(i+1);
            Tp00=Xr(i);
            Tp01=Xr(i-1);


            l = lambda(Tp00,ri);
            l2 = lambda2(Tp00,ri);

            a1=l;
            a2=l;
            a3=l2;

            Ar(i,i+1)   =   Ar(i,i+1)+ (a1/(dr^2)) + a2/(2*ri*dr);
            Ar(i,i-1)   =   Ar(i,i-1)+(a1/(dr^2)) - a2/(2*ri*dr);
            Ar(i,i)     =   Ar(i,i)+ -2*a1/(dr^2);
            br(i)       =   br(i) -a3*(Tp10^2-Tp01^2)/dr;
        end
        Ar(1,1)=1;
        Ar(1,2)=-1;
        l = lambda(Xr(end),r(end));
        Ar(nr+1,nr+1) = -(l/dr) - hconv;
        Ar(nr+1,nr) = l/dr;
        Xr2 = Ar\br;
        Xr0=Xr;
        Xr=Xr2;
    end    

    %% Volúmen de control V2
    C = 2*R_clad_ext / (R_agua^2 - R_clad_ext^2);
    
    for it=1:10
        Tam = mean(Xz(1:nz+1));
        
        if Tam>=Xr(end)
            Tam = min(Xz(1:nz+1));
        end
        Az=sparse(2*nz+2,2*nz+2);
        bz = zeros(2*nz+2,1);
        for i=1:nz+1
            if i==1

                Az(i,i) =1 ;
                Az(nz+1+i,nz+1+i) =1 ;
                bz(i)=T_inlet;
                bz(nz+1+i)=w0;
            
            
            else
                Tp00=Xz(i);
                wp=Xz(nz+1+i);

                b1=drho(Tp00);
                b2=rho(Tp00);
                c1=drho(Tp00)*h(Tp00) + dh(Tp00)*rho(Tp00);
                c2=h(Tp00)*rho(Tp00);


                Az(i,i)         =   Az(i,i)         +   wp*b1/dz ;
                Az(i,i-1)       =   Az(i,i-1)       -   wp*b1/dz ;
                Az(i,nz+1+i)    =   Az(i,nz+1+i)    +   b2/dz;
                Az(i,nz+1+i-1)  =   Az(i,nz+1+i-1)  -   b2/dz;

                Az(nz+1+i,i)        =   Az(nz+1+i,i)        +   wp*c1/dz  ;
                Az(nz+1+i,i-1)      =   Az(nz+1+i,i-1)      -   wp*c1/dz;
                Az(nz+1+i,nz+1+i)   =   Az(nz+1+i,nz+1+i)   +   c2/dz ;
                Az(nz+1+i,nz+1+i-1) =   Az(nz+1+i,nz+1+i-1) -   c2/dz;

                bz(nz+1+i) = bz(nz+1+i)+C*hconv*(Xr(end)-Tam);
            end

        end
        
        Xz2 = Az\bz;
        Xz0=Xz;
        Xz=Xz2;

    end
  globalerr=[globalerr,sum((Xr0-Xr).^2/(1+nr))];
  globalerrz1=[ globalerrz1,sum((Xz0(1:nz+1)-Xz(1:nz+1)).^2/(1+nz))];  
  globalerrz2=[ globalerrz2,sum((Xz0(nz+2:2*nz+2)-Xz(nz+2:2*nz+2)).^2/(1+nz))];
end

%% Gráficos

% figure
% 
% plot(r*1000,Xr,'b','linewidth',1.2)
% title('Perfil radial temperatura del núcleo')
% xlabel('r (mm)')
% ylabel('Ta (K)')
% grid on
% 
% figure
% plot(z,Xz(1:nz+1),'b','linewidth',1.2)
% title('Perfil axial temperatura del agua')
% xlabel('z (m)')
% ylabel('Ta (K)')
% grid on
% 
% figure
% plot(z,Xz(nz+1+1:end),'b','linewidth',1.2)
% title('Velocidad del agua')
% xlabel('z (m)')
% ylabel('v (m/s)')
% grid on
% 
% figure
% plot(globalerr,'b','linewidth',1.2)
% title('Error cuadrático por iteración de temperatura del núcleo')
% xlabel('Iteración')
% ylabel('Error cuadrático')
% grid on
% 
% figure
% plot(globalerrz1,'b','linewidth',1.2)
% title('Error cuadrático por iteración de temperatura del agua')
% xlabel('Iteración')
% ylabel('Error cuadrático')
% grid on
% 
% figure
% plot(globalerrz2,'b','linewidth',1.2)
% title('Error cuadrático por iteración de velocidad del agua')
% xlabel('Iteración')
% ylabel('Error cuadrático')
% grid on



tdoppler = 0.3*Xr(1) + 0.7*Xr(end);
T_final = mean(Xz(1:nz+1));
rho_final = rho(T_final);
%fprintf('La temperatura de Doppler es T_D = %s\n',tdoppler)
%fprintf('La temperatura del agua media es T_a = %s\n',T_final)
%fprintf('La densidad del agua media es rho_a = %s\n',rho_final)

%% Funciones

function [l] = lambda(T,r)
    R_UO2 = 8.239*0.001/2; %m
    e=0.571*0.001; %m
    R_clad_ext = 9.517*0.001/2; %m
    R_clad_int = R_clad_ext - e; %m
    if r < R_UO2
        l = 1.05 + 2150/(T-73.15);
    elseif (R_UO2 < r)&& (r< R_clad_int)
        l=e*10^4;
    else
        l= 7.51 + 2.09e-2*T -1.45e-5*T^2 + 7.67e-9*T^3;
    end
        
end

function [l] = lambda2(T,r)
    R_UO2 = 8.239*0.001/2; %m
    e=0.571*0.001; %m
    R_clad_ext = 9.517*0.001/2; %m
    R_clad_int = R_clad_ext - e; %m

    if r < R_UO2
        l = -2150/(T-73.15)^2;
    elseif (R_UO2 < r)&& (r< R_clad_int)
        l=0;
    else
        l= 2.09e-2 -2*1.45e-5*T + 3*7.67e-9*T^2;
    end
        
end

function [r] = rho(T)
    r =  (-1.946*10^(-2)*T^2 + 20.317*T - 4523.7);       
end

function [hh] = h(T)
    hh =  (2.341*10^(-2)*T^2 -21.455*T + 5944.5)*1000;     
end

function [r] = drho(T)
    r =  (-1.946*10^(-2)*T*2 + 20.317);       
end

function [hh] = dh(T)
    hh =  (2.341*10^(-2)*T*2 -21.455)*1000;     
end
