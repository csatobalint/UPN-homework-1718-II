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
alfa=sqrt(s/m);
D=k/(2*alfa*m);
warr=linspace(0,20,20000);
lambda=warr/alfa;
N=(sqrt((1-lambda.^2).^2+4*D^2*lambda.^2)).^(-1);
figure
plot(lambda,N);
ylim([0 100]) 



