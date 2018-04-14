function [U_next_L, F_next_L, v_next_L, p_next_L, rho_next_L, T_next_L] = Infl_res_LHS (U_prev, A, c_v, gamma, R, dx, dt, p_res)

% Unpacking the primitive variables at the previous unit time-step points:
    v_prev=U_prev(2,:)./U_prev(1,:);
    T_prev=(U_prev(3,:)./U_prev(1,:))/c_v -(v_prev.^2)/(2*c_v);
    p_prev=R*((U_prev(1,:)/A).*T_prev);
    a_prev=sqrt(gamma*R*T_prev);             % a & v values at the previous time step, in each x pts
    
% Compressible method of characteristics:

% Origin of the characteristic line
    x_R  = dt*(a_prev(1)-v_prev(1))   /   ...
            (1 + (v_prev(2)-v_prev(1))/dx*dt - (a_prev(2)-a_prev(1))/dx*dt );  
% Interpolating all important variables to x_R and X_L:
    v_R=v_prev(1)+(v_prev(2)-v_prev(1))/dx*x_R;
    a_R=a_prev(1)+(a_prev(2)-a_prev(1))/dx*x_R;
    
% The inlet boundary condition for the LHS:
    theta=1;
    p_ref=10^5;
    T_ref=293;
    %p_res=1.1*10^5;                 % res : reservoir
    T_res=T_ref*(p_res/p_ref)^((gamma-1)/gamma);
    a_t_res=sqrt(gamma*R*T_res);    
    % rho_res=p_res/(R*T_res);
    % c_p=1007;   % !!
    % p_tank=1.2*p_ref;
    % p_pipe=p_ref;
% The conserved quantity is:
    K=a_R-theta*(gamma-1)/2*v_R;
    AA=((gamma-1)/2)^2+(gamma-1)/2;
    BB=theta*(gamma-1)*K;
    CC=K^2-a_t_res^2;
    v_next_L=(-BB+theta*sqrt(BB^2-4*AA*CC))/(2*AA);
    msg = 'Warning: Inlet boundary condition is not fulfilled.';

% if v_next_L > 0
%if K < a_t_res
    % The new isentropic wave propagation speed:
    a_next_L=K+theta*(gamma-1)/2*v_next_L;
    T_next_L=a_next_L^2/(gamma*R);
    p_next_L=p_ref*(T_next_L/T_ref)^(gamma/(gamma-1));
    rho_next_L=p_next_L/(R*T_next_L);
    e_next_L=c_v*T_next_L+0.5*v_next_L^2;
%else
    %error(msg)
%end

%Updating the U_next & F_next matrices on the boundaries:
    U_next_L=zeros(3,1);
    U_next_L(1)=A*rho_next_L;
    U_next_L(2)=A*rho_next_L*v_next_L;     
    U_next_L(3)=A*rho_next_L*e_next_L;
       
    F_next_L=zeros(3,1);
    F_next_L(1)=rho_next_L*v_next_L*A;
    F_next_L(2)=(rho_next_L*(v_next_L)^2+p_next_L)*A;
    F_next_L(3)=(rho_next_L*v_next_L*e_next_L+p_next_L*v_next_L)*A;

end