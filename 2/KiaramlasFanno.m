function M_cso=KiaramlasFanno(kappa, Tt, pt, p0, L, D, lambda, R,dx)
% close all
% clc

%% FELADAT
% Egy tartályban szén-doxidot (kappa = 1.3) tárolunk Tt = 20 °C-os hõmérsékleten, pt = 1.7
% bar nyomáson. A tartályhoz L = 5m hosszú, D = 10mm belsõ átmérõjû acélcsõ kapcsolódik,
% mely hidraulikailag simának és tökéletesen hõszigeteltnek tekinthetõ. A csõ végén légköri nyomás
% uralkodik. Számítsa ki a tartályból távozó tömegáramot, valamint a közeg állapotjelzõit
% (nyomás, sûrûség, hõmérséklet, sebesség) a csõ elején és végén!

% Hiányzó adatok, amiket fel kell venni: gázállandó
%                                        csõsúrlódási tényezõ   

%% ADATOK
tic
%kappa = 1.4183;
% Tt = 20;   % °C
% pt = 1.3;  % bar
% p0 = 1;    % bar
% L = 4;     % m
% D = 50;    % mm
% 
% lambda = 0.02;
% R = 8314/29;                  % universal gas constant for dry air; % J/kg/K

M10 = 0.1; % kezdeti érték a Mach-számhoz

%% ÁTVÁLTÁSOK
Tt = Tt + 273;
pt = pt*1e5;
p0 = p0*1e5;
D = D*1e-3;

%% SZÁMÍTÁSOK

% Keresztmetszet
A = D^2*pi/4;

% A belépõ Mach-szám
options = optimset('Display','off');
M1 = fsolve(@(M1)machSzamKereses(M1,kappa,pt,p0,lambda,L,D),M10,options);
[M2, xmax] = machSzamCsoVegen(M1,kappa,lambda,L,D);
Fanno_param=zeros(1,4/dx+1);
M_cso=zeros(1,4/dx+1);
T_cso=zeros(1,4/dx+1);
P_cso=zeros(1,4/dx+1);
rho_cso=zeros(1,4/dx+1);
velocity_cso=zeros(1,4/dx+1);
P0_cso=zeros(1,4/dx+1);
fanno_cso=zeros(1,4/dx+1);
for kszi=1:(1+4/dx)
    Fanno_param(kszi)=lambda*(xmax-(kszi-1)*L/(4/dx))/D;
    [M_cso(kszi), T_cso(kszi), P_cso(kszi), rho_cso(kszi), velocity_cso(kszi), P0_cso(kszi), fanno_cso(kszi)] = flowfanno(kappa, Fanno_param(kszi), 'fannosub');
end

%Csõbéli Mach-szám alakulása

csobeosztas=0:dx:4;
%figure(1)
%plot(csobeosztas,M_cso)
%grid on
% Állapotjelzõk a csõ elején
T1 = Tt/(1 + (kappa-1)/2*M1^2);
p1 = (1 + (kappa-1)/2*M1^2)^(kappa/(1-kappa))*pt;
rho1 = p1/R/T1;
a1 = sqrt(kappa*R*T1);
v1 = M1*a1;
mPont1 = rho1*v1*A;

% Állapotjelzõk a csõ végén
T2 = Tt/(1 + (kappa-1)/2*M2^2);
p2 = p0;
rho2 = p2/R/T2;
a2 = sqrt(kappa*R*T2);
v2 = M2*a2;
mPont2 = rho2*v2*A;

% Az eredmények kiírása
fprintf('Adatok a csõ elején:\n')
fprintf('Mach-szám:    %.4f\n',M1)
fprintf('Nyomás:       %.4f bar\n',p1*1e-5)
fprintf('Hõmérséklet:  %.4f K\n',T1)
fprintf('Sûrûség:      %.4f kg/m^3\n',rho1)
fprintf('Sebesség:     %.4f m/s\n',v1)
fprintf('Hangsebesség: %.4f m/s\n',a1)
fprintf('Tömegáram:    %.5f kg/s\n',mPont1)
fprintf('\n')
fprintf('Adatok a csõ végén:\n')
fprintf('Mach-szám:    %.4f\n',M2)
fprintf('Nyomás:       %.4f bar\n',p2*1e-5)
fprintf('Hõmérséklet:  %.4f K\n',T2)
fprintf('Sûrûség:      %.4f kg/m^3\n',rho2)
fprintf('Sebesség:     %.4f m/s\n',v2)
fprintf('Hangsebesség: %.4f m/s\n',a2)
fprintf('Tömegáram:    %.5f kg/s\n',mPont2)
toc
end