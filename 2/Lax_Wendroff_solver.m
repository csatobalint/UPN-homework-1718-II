function [U_next_LW, F_next_LW, v_next_LW, p_next_LW, rho_next_LW, T_next_LW] = Lax_Wendroff_solver(U_prev, F_prev, ~, x, A, c_v, gamma, R, dx, dt, sp_pts)

% Interpolating the U & F values at half x points
% The points we're intrested in:
    %xq=1.5:dx:sp_pts;
    %xq=x(1)+0.5:dx:sp_pts;
    xq=(x(1)+x(2))/2:dx:(x(end)+x(end-1))/2;
% transzponálni kell az U & F mátrixot, mert így tud csak interpolálni. utána meg visszatranszponálni
    U_prev_interp=(interp1(x,U_prev',xq))';
    % F_prev_interp=(interp1(x,F_prev',xq))';

% Moving a half time-step-->getting the U values there 
    U_half=zeros(3,sp_pts-1);
    for ii=1:sp_pts-1
        U_half(:,ii)=U_prev_interp(:,ii)+dt/(2*dx)*(F_prev(:,ii)-F_prev(:,ii+1));
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

% Moving a unit time-step and calculating U values over there:
    U_next=zeros(3,sp_pts);
    for jj=2:sp_pts-1      % indexet csiszolni, kimenet U hogy nézzen ki?
        U_next(:,jj)=U_prev(:,jj)+dt/dx*(F_half(:,jj-1)-F_half(:,jj));
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

% Implementing the boundary contitions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% #1 Closed pipe BC
% % Unpacking the primitive variables at the previous unit time-step points:
%     v_prev=U_prev(2,:)./U_prev(1,:);
%     T_prev=(U_prev(3,:)./U_prev(1,:))/c_v -(v_prev.^2)/(2*c_v);
%     p_prev=R*((U_prev(1,:)/A).*T_prev);
%     a_prev=sqrt(gamma*R*T_prev);             % a & v values at the previous time step, in each x pts
%     
% % Compressible method of characteristics:
% 
% % Origin of the characteristic line
%     x_R  = dt*(a_prev(1)-v_prev(1))   /   ...
%             (1 + (v_prev(2)-v_prev(1))/dx*dt - (a_prev(2)-a_prev(1))/dx*dt ); 
%     x_NL = dt*(a_prev(end)+v_prev(end))   /   ...
%             (1 - (v_prev(end-1)-v_prev(end))/dx*dt - (a_prev(end-1)-a_prev(end))/dx*dt ); 
%     
% % Interpolating all important variables to x_R and X_L:
%     v_R=v_prev(1)+(v_prev(2)-v_prev(1))/dx*x_R;
%     a_R=a_prev(1)+(a_prev(2)-a_prev(1))/dx*x_R;
%     p_R=p_prev(1)+(p_prev(2)-p_prev(1))/dx*x_R;
%     T_R=T_prev(1)+(T_prev(2)-T_prev(1))/dx*x_R;
%     v_L=v_prev(end)+(v_prev(end-1)-v_prev(end))/dx*x_NL;
%     a_L=a_prev(end)+(a_prev(end-1)-a_prev(end))/dx*x_NL;
%     p_L=p_prev(end)+(p_prev(end-1)-p_prev(end))/dx*x_NL;
%     T_L=T_prev(end)+(T_prev(end-1)-T_prev(end))/dx*x_NL;
%     
% % Calculating the LHS properties:
%     theta=1;
%     a_next =a_R-theta*(gamma-1)/2*v_R;                      %  ?_R = a_r-(?-1)/2*v_R = a = ? 
%     T_next(1)=a_next.^2/(gamma*R);
%     v_next(1)=0;                                            % wall boundary condition
%     e_next(1)=c_v*T_next(1)+0.5*(v_next(1)).^2;             % BIZTOS KELL A SEBESSÉGES TAG?
%     p_next(1)=p_prev(1)*(T_next(1)/T_prev(1))^((gamma)/gamma-1);        % isentropic process along ?=const line
%     rho_next(1)=p_next(1)/(R*T_next(1));
%     
% % Calculating the RHS properties:
%     theta=-1;
%     a_next =a_L-theta*(gamma-1)/2*v_L;                      %  ?_R = a_r-(?-1)/2*v_R = a = ? 
%     T_next(end)=a_next.^2/(gamma*R);
%     v_next(end)=0;                                            % wall boundary condition
%     e_next(end)=c_v*T_next(end)+0.5*(v_next(end)).^2;             
%     p_next(end)=p_prev(end)*(T_next(end)/T_prev(end))^((gamma)/gamma-1);        % isentropic process on the boundary in time
%     rho_next(end)=p_next(end)/(R*T_next(end));    
%     
% 
% %Updating the U_next & F_next matrices on the boundaries:
%     U_next(1,:)=A*rho_next;
%     U_next(2,:)=A*rho_next.*v_next;     
%     U_next(3,:)=A*rho_next.*e_next;
%     F_next(1,:)=rho_next.*v_next*A;
%     F_next(2,:)=(rho_next.*(v_next).^2+p_next)*A;
%     F_next(3,:)=(rho_next.*v_next.*e_next+p_next.*v_next)*A;
%     

% %% #1 Inow from the reservoir
% % LHS
% theta=1;
% p_ref=10^5;
% T_ref=293;
% p_res=1.2*10^5;                 % res : reservoir
% T_res=T_ref*(p_res/p_ref)^((gamma-1)/gamma);
% a_t_res=sqrt(gamma*R*T_res);    
% 
% rho_res=p_res/(R*T_res);
% c_p=1007;   % !!
% % p_tank=1.2*p_ref;
% % p_pipe=p_ref;
% 
% K=a_R-theta*(gamma-1)/2*v_R;
% A=((gamma-1)/2)^2+(gamma-1)/2;
% B=theta*(gamma-1)*K;
% C=K^2-a_t_res^2;
% v_next_L=(-B+theta*sqrt(B^2-4*A*C))/(2*A);
% if v_next_L > 0
% % The new isentropic wave propagation speed:
% a_next=K+theta*(gamma-1)/2*v_next;
% T_next=a_next/(gamma*R);
% p_next=p_ref*(T_next/T_ref)^(gamma/(gamma-1));
% rho_next=p_next/(R*T_next);
% e_next=c_v*T_next+0.5*v_next^2;
% msg = 'Warning: Inlet boundary condition is not fulfilled.'
% else error(msg)

end
