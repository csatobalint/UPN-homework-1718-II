function fvbrRK4(T,N,m,k,s,f0,w,c1,c2)

h=T/N; % stepsize
t=linspace(0,T,N+1); % grid
y=zeros(2,N+1); % numerical solution
epsilon=10^-4;

% RK4:
y(:,1)=[c1 c2]';
% corr=1;
sol4=zeros(2,N+1);
sol5=zeros(2,N+1);

for j=1:N
    %h = min(h, 2-t);
    Y_1      = h*fvbr(t(j),y(:,j),m,k,s,f0,w);
    Y_2      = h*fvbr(t(j)+1/4*h, y(:,j)+1/4*Y_1,m,k,s,f0,w);
    Y_3      = h*fvbr(t(j)+3/8*h, y(:,j)+3/32*Y_1+9/32*Y_2,m,k,s,f0,w);
    Y_4      = h*fvbr(t(j)+12/13*h, y(:,j)+1932/2197*Y_1-7200/2197*Y_2+7296/2197*Y_3           ,m,k,s,f0,w);
    Y_5      = h*fvbr(t(j)+h, y(:,j)+439/216*Y_1-8*Y_2+3680/513*Y_3-845/4104*Y_4               ,m,k,s,f0,w);
    Y_6      = h*fvbr(t(j)+1/2*h, y(:,j)+8/27*Y_1-2*Y_2+3544/2565*Y_3+1859/4104*Y_4-11/40*Y_5  ,m,k,s,f0,w);
    
    
    % 4th order approximation:
    sol4(:,j+1) = y(:,j) + 25/216*Y_1+1408/2565*Y_3+2197/4101*Y_4-1/5*Y_5;
    % 5th order approximation:
    sol5(:,j+1) = y(:,j) + 16/135*Y_1+6656/12.825*Y_3+28.561/56.430*Y_4-9/50*Y_5+2/55*Y_6;
    % Error in solution (for y(:,j+1) between 4th & 5th order RK.)
%     error = abs(max(sol5(:,j+1)-y(:,j+1)));
    R = abs(max(sol5(:,j+1)-y(:,j+1)))/h;   
    delta = 0.84*(epsilon/R)^(1/4);
    
    if R<=epsilon
        t = t+h;
        y(:,j+1) = sol4(:,j+1);
        j=j+1;
        h = delta*h;
    else 
        h = delta*h;
    end
    
%     if      10^(-4) < error < 10^(-3)    % The solution is proper, we can move on
%         corr=1;
%         
%     elseif  error > 10^(-3)          % The sol. is inaccurate, we're halving the time-step & recalc. the sol. Then move on
%         corr=0.85;
%         
%         Y_1      = h*fvbr(t(j),y(:,j)                                                              ,m,k,s,f0,w);
%         Y_2      = h*fvbr(t(j)+1/4*corr*h, y(:,j)+1/4*Y_1                                               ,m,k,s,f0,w);
%         Y_3      = h*fvbr(t(j)+3/8*corr*h, y(:,j)+3/32*Y_1+9/32*Y_2            ,m,k,s,f0,w);
%         Y_4      = h*fvbr(t(j)+12/13*corr*h, y(:,j)+1932/2197*Y_1-7200/2197*Y_2+7296/2197*Y_3           ,m,k,s,f0,w);
%         Y_5      = h*fvbr(t(j)+corr*h, y(:,j)+439/216*Y_1-8*Y_2+3680/513*Y_3-845/4104*Y_4               ,m,k,s,f0,w);
%         Y_6      = h*fvbr(t(j)+1/2*corr*h, y(:,j)+8/27*Y_1-2*Y_2+3544/2565*Y_3+1859/4104*Y_4-11/40*Y_5  ,m,k,s,f0,w);
% %       y(:,j+1) = y(:,j) + 25/216*Y_1+1408/2565*Y_3+2197/4101*Y_4-1/5*Y_5;
%         y(:,j+1) = y(:,j) + 16/135*Y_1+6656/12.825*Y_3+28.561/56.430*Y_4-9/50*Y_5+2/55*Y_6;
%                 
%     else                            % The sol. is too accurate, we're duplicating the time-step but no recalc.
%         corr=1.2;
%         
%     end
%         j=j+1;
       
    
end


figure
subplot(1,2,1)
plot(t,y(1,:),'red')
xlabel('Time')
subplot(1,2,2)
plot(y(1,:),y(2,:))
title('Phase portrait')
% y(1,end)

end