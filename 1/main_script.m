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

%initial conditions
c1=2;
c2=2;

%timestep
%h=0.005;
%N=(T-t0)/h;
N=20000;


%% forced vibration

%[amp, amp_data, time] = fvbrRK4(T,N,m,k,s,f0,w,c1,c2)

%% For loop for the amplification diagram

ii=1;
for w = 2:0.5:10
    amp(ii)= fvbrRK4_amp(T,N,m,k,s,f0,w,c1,c2);
    w_list(ii)=w;
    ii=ii+1;
end

figure
plot(w_list/sqrt(s/m), amp*s/f0, '*')
title('Amplification diagram')
xlabel('$\lambda=\frac{\omega}{\alpha}$','interpreter','latex')
xlabel('$N$','interpreter','latex')

%%
fvbrRK4(T,N,m,k,s,f0,w,c1,c2)