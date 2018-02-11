%% UNP-homework-I

clear variables;
close all;
clc;
format long

%time range
t0=0;
T=120;

%parameters of the mechanical system
m=2;
k=0.2;
s=100;
f0=60;
w=5;

%initial conditions
c1=2;
c2=2;

%timestep
h=0.01;
N=(T-t0)/h;

%forced vibration
fvbrRK4(T,N,m,k,s,f0,w,c1,c2)
