function [U_next_LW, F_next_LW, Q_next_LW, v_next_LW, p_next_LW, rho_next_LW, T_next_LW, e_next_LW] = Lax_Wendroff_solver(U_prev, F_prev, Q_prev, x, A, c_v, gamma, R, dx, dt, sp_pts, lambda, d)

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
    Q_half(2,:)=-A/2*lambda/d*rho_half.*v_half.*abs(v_half);

% Moving a unit time-step and calculating U values over there:
    U_next=zeros(3,sp_pts);
    for jj=2:sp_pts-1      % indexet csiszolni, kimenet U hogy nézzen ki?
        U_next(:,jj)=U_prev(:,jj)+dt/dx*(F_half(:,jj-1)-F_half(:,jj))+Q_half(:,jj);
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

end
