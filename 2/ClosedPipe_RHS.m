function [U_next_R, F_next_R, v_next_R, p_next_R, rho_next_R, T_next_R] = ClosedPipe_RHS (U_prev, A, c_v, gamma, R, dx, dt)

% Unpacking the primitive variables at the previous unit time-step points:
    v_prev=U_prev(2,:)./U_prev(1,:);
    T_prev=(U_prev(3,:)./U_prev(1,:))/c_v -(v_prev.^2)/(2*c_v);
    p_prev=R*((U_prev(1,:)/A).*T_prev);
    a_prev=sqrt(gamma*R*T_prev);             % a & v values at the previous time step, in each x pts
    
% Compressible method of characteristics:

% Origin of the characteristic line
    x_NL = dt*(a_prev(end)+v_prev(end))   /   ...
            (1 - (v_prev(end-1)-v_prev(end))/dx*dt - (a_prev(end-1)-a_prev(end))/dx*dt ); 
    
% Interpolating all important variables to x_R and X_L:
    v_L=v_prev(end)+(v_prev(end-1)-v_prev(end))/dx*x_NL;
    a_L=a_prev(end)+(a_prev(end-1)-a_prev(end))/dx*x_NL;
%     p_L=p_prev(end)+(p_prev(end-1)-p_prev(end))/dx*x_NL;
%     T_L=T_prev(end)+(T_prev(end-1)-T_prev(end))/dx*x_NL;
       
% Calculating the RHS properties:
    theta=-1;
    a_next_R =a_L-theta*(gamma-1)/2*v_L;                      %  ?_R = a_r-(?-1)/2*v_R = a = ? 
    T_next_R=a_next_R.^2/(gamma*R);
    v_next_R=0;                                            % wall boundary condition
    e_next_R=c_v*T_next_R+0.5*(v_next_R).^2;             
    p_next_R=p_prev(end)*(T_next_R/T_prev(end))^((gamma)/gamma-1);        % isentropic process on the boundary in time
    rho_next_R=p_next_R/(R*T_next_R);    
    
%Updating the U_next & F_next matrices on the boundaries:
    U_next_R=zeros(3,1);
    U_next_R(1)=A*rho_next_R;
    U_next_R(2)=A*rho_next_R*v_next_R;     
    U_next_R(3)=A*rho_next_R*e_next_R;
    
    F_next_R=zeros(3,1);
    F_next_R(1)=rho_next_R*v_next_R*A;
    F_next_R(2)=(rho_next_R*(v_next_R)^2+p_next_R)*A;
    F_next_R(3)=(rho_next_R*v_next_R*e_next_R+p_next_R*v_next_R)*A;

