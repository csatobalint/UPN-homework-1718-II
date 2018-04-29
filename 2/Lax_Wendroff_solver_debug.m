%% Applying Lax-Wendroff scheme (and other methods) for gas outflow from a tank.
% Unsteady flows in pipe networks home assignment 
% The script was made by Zsigmond Zalán and Csató Bálint.
%     clear all
%     close all
%     clc
    % format shortEng
    % format long
    format compact
    format short

% Material, geometrical properties
    R=8314/29;                  % universal gas constant for dry air
    d=0.05;                     % pipe diameter, [m]
    A=d^2*pi/4;                 % [m^2]
    c_v=710;                    % [J/(kg*K)] %[kJ/kg/K] ftp://ftp.energia.bme.hu/pub/Tuzelestechnika/MSc/fajho.pdf  
    c_p=1007;                   % [J/(kg*K)]
    gamma=c_p/c_v;
    lambda=0.02;
    
% Pipe geometry settings    
    xL=0;                       % left side of the pipe
    xR=4;                      % right side of the pipe
    L=xR-xL;                    % length of the pipe
    
% Spatial discretization settings 
    dx=0.2;                       % spatial step size
    %sp_pts=(L+1)/dx            % spatial resolution
    x=xL:dx:xR;                 % spatial grid
    sp_pts=length(x);
    xq=(x(1)+x(2))/2:dx:(x(end)+x(end-1))/2;
    
% Initialization of state variables  
    p_ref=10^5;
    p_res=1.1*10^5;
    p=p_ref*ones(1,sp_pts);      % [Pa]
    p(floor(end/2))=1.1*10^5;   % pressure peak at the middle
    T=300*ones(1,sp_pts);       % [K]
    v=zeros(1,sp_pts);          % velocity field [m/s]
    rho=p./(R*T);               % density (ideal gas law) [kg/m3]
    e=c_v*T+(v.^2)/2;           % internal energy [J]
    %estimation of variables
    v_est=sqrt((p_res-p_ref)*2/mean(rho));

% Time discretization settings 
    a=sqrt(gamma*R*T(1,1));
    CFL=0.6;
    dt=CFL*dx/(a+v(1));                    % dx/dt<=a --> dx/dt:=0.5*a --> dt~0.005
    tsteps=1000;
    t=zeros(1,tsteps);
    
% Initialization of state vectors
    % Creating empty arrays
    U_ini=zeros(3,sp_pts);
    F_ini=zeros(3,sp_pts);
    Q_ini=zeros(3,sp_pts);
    
    %Calculating starting values
    U_ini(1,:)=A*rho;
    U_ini(2,:)=A*rho.*v;
    U_ini(3,:)=A*rho.*e;
    
    F_ini(1,:)=A*rho.*v;
    F_ini(2,:)=A*(rho.*v.^2+p);
    F_ini(3,:)=A*(rho.*e.*v+p.*v);
    
%   F_s=A*rho/2*lambda/D.*v.*abs(v);
%     F_s_ini=0;
        
    U_prev=U_ini;
    F_prev=F_ini;
    Q_prev=Q_ini;
    
    % Storage of U & F matrices:
    U_write=zeros(3,sp_pts,tsteps);
    F_write=zeros(3,sp_pts,tsteps);
    v_write=zeros(tsteps,sp_pts);
    p_write=zeros(tsteps,sp_pts);
    rho_write=zeros(tsteps,sp_pts);
    T_write=zeros(tsteps,sp_pts);
    t_write=zeros(tsteps,sp_pts);
    dv_write=zeros(1,tsteps);

% Interpolating the U & F values at half x points
% The points we're intrested in:
    %xq=1.5:dx:sp_pts;
    %xq=x(1)+0.5:dx:sp_pts;
    xq=(x(1)+x(2))/2:dx:(x(end)+x(end-1))/2;
% transzponálni kell az U & Q mátrixot, mert így tud csak interpolálni. utána meg visszatranszponálni
    U_prev_interp=(interp1(x,U_prev',xq))';
    Q_prev_interp=(interp1(x,Q_prev',xq))';

% Moving a half time-step-->getting the U values there 
    U_half=zeros(3,sp_pts-1);
    for ii=1:sp_pts-1
        U_half(:,ii)=U_prev_interp(:,ii)+dt/(2*dx)*(F_prev(:,ii)-F_prev(:,ii+1))+dt/2*Q_prev_interp(:,ii);
    end

% Unpacking the primitive variables at half time-step points:
    rho_half=U_half(1,:)/A;
    v_half=U_half(2,:)./U_half(1,:);
    e_half=U_half(3,:)./U_half(1,:);
    T_half=1/c_v*(e_half-(v_half.^2)/2);
    p_half=R*(rho_half.*T_half);
    
% Calculating the F values at half time-step points:
    F_half(1,:)= rho_half.*v_half*A;
    F_half(2,:)=(rho_half.*(v_half).^2+p_half)*A;
    F_half(3,:)=(rho_half.*v_half.*e_half+p_half.*v_half)*A;
    Q_half=zeros(3,sp_pts-1);
    Q_half(2,:)=-A*rho_half./2*lambda/d.*v_half.*abs(v_half);

% Moving a unit time-step and calculating U values over there:
    U_next=zeros(3,sp_pts);
    for ii=2:sp_pts-1      % indexet csiszolni, kimenet U hogy nézzen ki?
        U_next(:,ii)=U_prev(:,ii)+dt/dx*(F_half(:,ii-1)-F_half(:,ii))+Q_half(:,ii);
    end
% The solver doesn't work at the boundaries so we can cut out the bou. pts:
    U_next_LW=U_next(:,2:end-1);
% Unpacking the primitive variables at the next unit time-step points:
    rho_next_LW=U_next_LW(1,:)/A;
    v_next_LW=U_next_LW(2,:)./U_next_LW(1,:);
    e_next_LW=U_next_LW(3,:)./U_next_LW(1,:);
    T_next_LW=1/c_v*(e_next_LW-(v_next_LW.^2)./2);
    p_next_LW=R*(rho_next_LW.*T_next_LW);
    
% Calculating the F values at the next unit time-step points:
    F_next_LW(1,:)=rho_next_LW.*v_next_LW*A;
    F_next_LW(2,:)=(rho_next_LW.*(v_next_LW).^2+p_next_LW)*A;
    F_next_LW(3,:)=(rho_next_LW.*v_next_LW.*e_next_LW+p_next_LW.*v_next_LW)*A;
    Q_next_LW=zeros(3,sp_pts-2);
    Q_next_LW(2,:)=-A/2*lambda/d*rho_next_LW.*v_next_LW.*abs(v_next_LW);
