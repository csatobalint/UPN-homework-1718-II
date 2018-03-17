%% Unsteady flows in pipe networks homework. Applying Lax-Wendroff scheme (and other methods) for gas outflow from a tank.
clear all
close all
clc
% format shortEng
% format long
format compact
format short

R=8314/29;                  % universal gas constant for dry air
d=0.05;
A=d^2*pi/4;                        % [m^2]
c_v=710;                    % [J/(kg*K)] %[kJ/kg/K] ftp://ftp.energia.bme.hu/pub/Tuzelestechnika/MSc/fajho.pdf  
c_p=1007;                   % [J/(kg*K)]
gamma=c_p/c_v;
sp_pts=50;
x=1:sp_pts;
% p=10^5*[10:-0.5:5.5];
p=10^5*ones(1,sp_pts);      % [Pa]
% p(20)=2*p(20);
% p(21)=2*p(21);
p(floor(end/2))=1.5*10^5;    
T=300*ones(1,sp_pts);       % [K]
v=zeros(1,sp_pts);
rho=p./(R*T);               % ideal gas law
e=c_v*T+(v.^2)/2;

dx=1;
a=sqrt(gamma*R*T(1,1));
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

U_prev=U_ini;
F_prev=F_ini;
Q_prev=Q_ini;
tsteps=5;
% storage of U & F matrices:
U_write=zeros(3,sp_pts,tsteps);
F_write=zeros(3,sp_pts,tsteps);
v_write=zeros(tsteps,sp_pts);
p_write=zeros(tsteps,sp_pts);
for i=1:tsteps  % moving in the time domain
    [U_next, F_next,v_next, p_next] = Lax_Wendroff_solver(U_prev, F_prev, Q_prev, x, A, c_v, gamma, R, dx, dt, sp_pts);
%     U_next(isnan(U_next))=0;
%     F_next(isnan(F_next))=0;
    U_write(:,:,i)=U_next;
    F_write(:,:,i)=F_next;
    v_write(i,:)=v_next;
    p_write(i,:)=p_next;
    U_prev=U_next; F_prev=F_next;
end
 
% figure(3)
% p_write(isnan(p_write))=0;
% plot(x,p,'-*k', x,p_write(1,:),'-or')
% legend('initial step','1st timestep','2nd timestep','Location','northeast')

figure(1)
v_write(isnan(v_write))=0;
plot(x,v,'-*k', x,v_write(1,:),'-or',x, v_write(2,:),'-*g', x, v_write(3,:),'-db', x, v_write(4,:),'-+c', x, v_write(5,:),'-ok')
legend('initial step','1st timestep','2nd timestep','3rd timestep','4th timestep','5th timestep','Location','northeast')
title('Velocity')
figure(2)
p_write(isnan(p_write))=0;
plot(x,p,'-*k', x,p_write(1,:),'-or',x, p_write(2,:),'-*g', x, p_write(3,:),'-db', x, p_write(4,:),'-+c', x, p_write(5,:),'-ok')
legend('initial step','1st timestep','2nd timestep','3rd timestep','4th timestep','5th timestep','Location','northeast')
title('Pressure')

U_write
F_write