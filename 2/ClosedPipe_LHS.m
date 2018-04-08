function [U_next_L, F_next_L, v_next_L, p_next_L, rho_next_L, T_next_L] = ClosedPipe_LHS (U_prev, A, c_v, gamma, R, dx, dt)

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
%     p_R=p_prev(1)+(p_prev(2)-p_prev(1))/dx*x_R;
%     T_R=T_prev(1)+(T_prev(2)-T_prev(1))/dx*x_R;

% Calculating the LHS properties:
    theta=1;
    a_next_L =a_R-theta*(gamma-1)/2*v_R;                      %  ?_R = a_r-(?-1)/2*v_R = a = ? 
    T_next_L=a_next_L.^2/(gamma*R);
    v_next_L=0;                                            % wall boundary condition
    e_next_L=c_v*T_next_L+0.5*(v_next_L).^2;             
    p_next_L=p_prev(1)*(T_next_L/T_prev(1))^((gamma)/gamma-1);        % isentropic process along ?=const line
    rho_next_L=p_next_L/(R*T_next_L);   

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