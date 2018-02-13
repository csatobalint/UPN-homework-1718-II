%% UNP-homework-I

%%
clear all;
close all;
clc;
format long

%time range
t0=0;
T=400;

%parameters of the mechanical system
m=2;
k=0.2;
s=100;
f0=60;
w=10;
warr=linspace(0,20,1000);
lambda=warr/sqrt(s/m);
D=k/(2*sqrt(s/m)*m);
N=(sqrt((1-lambda.^2).^2+4*D^2*lambda.^2)).^(-1);
plot(lambda,N);


%%
%initial conditions
c1=2;
c2=2;

%timestep
%h=0.005;
%N=(T-t0)/h;
N=5000;

%% forced vibration

amp = fvbrRK4(T,N,m,k,s,f0,w,c1,c2)

%% For loop for the amplification diagram

ii=1;
for w = 2:0.5:10
    amp(ii)= fvbrRK4(T,N,m,k,s,f0,w,c1,c2);
    w_list(ii)=w;
    ii=ii+1;
end

figure 
plot(w_list, amp, '*')

