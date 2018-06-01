%% Applying Lax-Wendroff scheme (and other methods) for gas outflow from a tank.
% Unsteady flows in pipe networks home assignment 
% The script was made by Zsigmond Zalán and Csató Bálint.
    clear all
    close all
    clc
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
    dx=0.01;                       % spatial step size
    %sp_pts=(L+1)/dx            % spatial resolution
    x=xL:dx:xR;                 % spatial grid
    sp_pts=length(x);
    xq=(x(1)+x(2))/2:dx:(x(end)+x(end-1))/2;
    
% Initialization of state variables  
    p_ref=10^5;
    p_res=1.3*10^5;
    p=p_ref*ones(1,sp_pts);      % [Pa]
    %p(floor(end/2))=1.1*10^5;   % pressure peak at the middle
%     p([1:floor(end/2)])=1.1*10^5;
    period = 2*L;
    %p=p+0.1*10^5*sin(2*3.14/period*x);
    T=300*ones(1,sp_pts);       % [K]
    v=zeros(1,sp_pts);          % velocity field [m/s]
    rho=p./(R*T);               % density (ideal gas law) [kg/m3]
    e=c_v*T+(v.^2)/2;           % internal energy [J]
    %estimation of variables
    v_est=sqrt((p_res-p_ref)*2/mean(rho));

% Time discretization settings 
    a=sqrt(gamma*R*T(1,1));
    CFL=0.8;
    dt=CFL*dx/(a+v(1));                    % dx/dt<=a --> dx/dt:=0.5*a --> dt~0.005
    tsteps=10000;
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
tic       
last_timestep=tsteps;
for i=1:tsteps  % moving in the time domain
    
    [U_next_LW, F_next_LW, Q_next_LW, v_next_LW, p_next_LW,rho_next_LW, T_next_LW, e_next_LW] = Lax_Wendroff_solver(U_prev, F_prev, Q_prev, x, A, c_v, gamma, R, dx, dt, sp_pts, lambda, d);
    [U_next_L, F_next_L, v_next_L, p_next_L, rho_next_L, T_next_L] = Infl_res_LHS(U_prev, A, c_v, gamma, R, dx, dt, p_res);
    [U_next_R, F_next_R, v_next_R, p_next_R, rho_next_R, T_next_R] = Outfl_RHS(U_prev, A, c_v, gamma, R, dx, dt, p_ref);
%     [U_next_L, F_next_L, v_next_L, p_next_L, rho_next_L, T_next_L] = ClosedPipe_LHS(U_prev, A, c_v, gamma, R, dx, dt);
%     [U_next_R, F_next_R, v_next_R, p_next_R, rho_next_R, T_next_R] = ClosedPipe_RHS(U_prev, A, c_v, gamma, R, dx, dt);
    
    
    U_next=[U_next_L U_next_LW U_next_R];
    F_next=[F_next_L F_next_LW F_next_R];
    v_next=[v_next_L v_next_LW v_next_R];
    p_next=[p_next_L p_next_LW p_next_R];
    rho_next=[rho_next_L rho_next_LW rho_next_R];
    T_next=[T_next_L T_next_LW T_next_R];
    %     U_next(isnan(U_next))=0;
    %     F_next(isnan(F_next))=0;
    
%     %OnlyWalls
%      [U_next, F_next, Q_next, v_next, p_next, rho_next, T_next] = Lax_Wendroff_solver_wall(U_prev, F_prev, Q_prev, x, A, c_v, gamma, R, dx, dt, sp_pts, lambda, d);

    
    % Write the results of the current time step
        U_write(:,:,i)=U_next;
        F_write(:,:,i)=F_next;
        v_write(i,:)=v_next;
        p_write(i,:)=p_next;
        rho_write(i,:)=rho_next;
        T_write(i,:)=T_next;
        if i>1
        t(i)=t(i-1)+dt;
        else
        t(i)=dt;
        end
        
     % Calculate residuals
        v_prev=U_prev(2,:)./U_prev(1,:);
        dv_next=max(   abs(   (v_next-v_prev)/max(max(v_prev),max(v_next)))    );
        if dv_next<1*10^-5
        %if dv_next<1*10^-2    %wall!!!!!! 
            last_timestep=i;
            break 
        end
        dv_write(i)=dv_next;
        
    % Initialize the next time step
        U_prev=U_next; 
        F_prev=F_next;
        
        t_write=CFL*dx./(v_next+sqrt(gamma*R*T_next));
        dt=min(CFL*dx./(v_next+sqrt(gamma*R*T_next)));
    
    fprintf('Step %d: t = %6.4f, dt = %10.7f\n',i, t(i),dt);
    
    p_final=U_write(2,:,tsteps)./U_write(1,:,tsteps);
    
end
toc

semilogy(t,dv_write);
xlabel('$t~[s]$','interpreter','latex');
ylabel('$\delta v~[-]$','interpreter','latex');
title('Relative changes of velocity');
grid on
set(gca,'FontSize',14)
filename = ['plots/residuals_dx_' num2str(dx) '.png'];
%filename = ['plots/residuals_wall_' num2str(dx) '.png'];
saveas(gcf,filename)
close all
tic