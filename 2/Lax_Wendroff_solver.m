function [U_next, F_next] = Lax_Wendroff_solver(U_ini, F_ini, Q_ini, x, A, c_v, R, dx, dt, p, T, v, rho)

% Interpolating the U & F values at half x points
% The points we're intrested in:
xq=1.5:10;
% transzponálni kell az U & F mátrixot, mert így tud csak interpolálni. utána meg visszatranszponálni
U_prev_interp=(interp1(x,U_ini',xq))'; 
% F_prev_interp=(interp1(x,F_ini',xq))';

% Moving a half time-step-->getting the U values there //delta(t)-t és
% delta(x)-et majd deklarálni kell
U_half=zeros(3,9);
for ii=1:9
    U_half(:,ii)=U_prev_interp(:,ii)+dt/(2*dx)*(F_ini(:,ii)-F_ini(:,ii+1));
end

% Calculating the F values at half time-step points:
rho_half=U_half(1,:)/A;
v_half=U_half(2,:)/U_half(1,:);
e_half=U_half(3,:)/U_half(1,:);
T_half=1/c_v*(e_half-(v_half.^2)/2);
p_half=R*(rho_half.*T_half);
F_half(1,:)=rho_half.*v_half*A;
F_half(2,:)=(rho_half.*(v_half).^2+p_half)*A;
F_half(3,:)=(rho_half.*v_half.*e_half+p_half.*v_half)*A;

% Moving a unit time-step:
U_next=zeros(3,10);
for jj=2:9      % indexet csiszolni, kimenet U hogy nézzen ki?
    U_next(:,jj)=U_ini(:,jj)+dt/dx*(F_half(:,jj-1)-F_half(:,jj));
end

% Calculating the F values at the unit time-step points:
rho_next=U_next(1,:)/A;
v_next=U_next(2,:)./U_next(1,:);
e_next=U_next(3,:)./U_next(1,:);
T_next=1/c_v*(e_next-(v_next.^2)./2);
p_next=R*(rho_next.*T_next);
F_next(1,:)=rho_next.*v_next*A;
F_next(2,:)=(rho_next.*(v_next).^2+p_next)*A;
F_next(3,:)=(rho_next.*v_next.*e_next+p_next.*v_next)*A;

%%
subplot(2,2,1)
plot(x, p_next, x, p, 'r')
subplot(2,2,2)
plot(x, v_next, x, v, 'r')
subplot(2,2,3)
plot(x, rho_next, x, rho, 'r')
subplot(2,2,4)
plot(x, T_next, x, T, 'r')


end
