%% UNP-homework-I

clear all;
close all;
clc;
format long

%time range
t0=0;
T=200;

%parameters of the mechanical system
m=2;
k=0.2;
s=100;
f0=60;
w=2;
F0=f0/s;

%initial conditions
c1=2;
c2=2;

N=20000;
figure(1)
amp = fvbrRK4(T,N,m,k,s,f0,w,c1,c2)
%% For loop for the amplification diagram

ii=1;
for w = 0:0.05:14
    amp(ii)= fvbrRK4_amp(T,N,m,k,s,f0,w,c1,c2);
    w_list(ii)=w;
    ii=ii+1;
end

alfa=sqrt(s/m);
D=k/(2*alfa*m);
warr=linspace(0,14,10000);
lambda=warr/alfa;
N=(sqrt((1-lambda.^2).^2+4*D^2*lambda.^2)).^(-1);

% Numerical solution
num_sol_x=w_list/sqrt(s/m);
num_sol_y=amp*s/f0;
figure(2)
semilogy(num_sol_x(2:end), num_sol_y(2:end), '*')
% plot(num_sol_x(2:end), num_sol_y(2:end), '*')
title('Amplification diagram with numerical solution')
xlabel('$\lambda=\frac{\omega}{\alpha}$','interpreter','latex')
ylabel('$N$','interpreter','latex')
ylim([0 100])

% Analytical solution
figure(3)
semilogy(lambda,N);
% plot(lambda,N);
ylim([0 100])
title('Amplification diagram with analytical solution')
xlabel('$\lambda=\frac{\omega}{\alpha}$','interpreter','latex')
ylabel('$N$','interpreter','latex')
ylim([0 100])
