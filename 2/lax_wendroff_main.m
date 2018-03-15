%% Unsteady flows in pipe networks homework. Applying Lax-Wendroff scheme (and other methods) for gas outflow from a tank.
clear all
close all
clc
format shortEng
% format long
format compact

R=8314/29;                  % universal gas constant for dry air
d=0.05;
A=d^2*pi/4;                        % [m^2]
c_v=710;                    % [J/(kg*K)] %[kJ/kg/K] ftp://ftp.energia.bme.hu/pub/Tuzelestechnika/MSc/fajho.pdf  
c_p=1007;                   % [J/(kg*K)]
sp_pts=10;
x=1:sp_pts;
% p=10^5*[10:-0.5:5.5];
p=10^5*ones(1,sp_pts);      % [Pa]
p(5)=2*p(5);
p(6)=2*p(6);
% p(1:floor(end/2))=1.1*10^5;      %
T=300*ones(1,sp_pts);       % [K]
v=zeros(1,sp_pts);
rho=p./(R*T);                 % ideal gas law
e=c_v*T+(v.^2)/2;

dx=1;
a=sqrt(c_p/c_v*R*T(1,1));
dt=dx/a;                  % dx/dt<=a --> dx/dt:=0.5*a --> dt~0.005

U_ini=zeros(3,sp_pts);
F_ini=zeros(3,sp_pts);
Q_ini=zeros(1,sp_pts);

U_ini(1,:)=A*rho;
U_ini(2,:)=A*rho.*v;
U_ini(3,:)=A*rho.*e;
F_ini(1,:)=A*rho.*v;
F_ini(2,:)=A*(rho.*v.^2+p);
F_ini(3,:)=A*(rho.*e.*v+p.*v);

[U_next, F_next] = Lax_Wendroff_solver(U_ini, F_ini, Q_ini, x, A, c_v, R, dx, dt, p, T, v, rho);

% for i=1:20  % moving in the time domain
%     % itt kéne meghívni a függvényemet!
% end
U_ini
U_next
F_ini
F_next
    

%% nterp ellenõrzésére
xq=1.5:10;
U_prev_interp=(interp1(x,U_ini',xq))';
U_ini
U_prev_interp