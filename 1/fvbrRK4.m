function fvbrRK4(T,N,m,k,s,f0,w,c1,c2)

h=T/N; % stepsize
% t=linspace(0,T,N+1); % grid
y=zeros(2,N+1); % numerical solution
epsilon=10^-5;

% RK4:
y(:,1)=[c1 c2]';
t=0;
% corr=1;
sol4=zeros(2,N+1);
sol5=zeros(2,N+1);

for ii=1:N
    
    Y_1      = h*fvbr(t,y(:,ii)                                                              ,m,k,s,f0,w);
    Y_2      = h*fvbr(t+1/4*h, y(:,ii)+1/4*Y_1                                               ,m,k,s,f0,w);
    Y_3      = h*fvbr(t+3/8*h, y(:,ii)+3/32*Y_1+9/32*Y_2                                     ,m,k,s,f0,w);
    Y_4      = h*fvbr(t+12/13*h, y(:,ii)+1932/2197*Y_1-7200/2197*Y_2+7296/2197*Y_3           ,m,k,s,f0,w);
    Y_5      = h*fvbr(t+h, y(:,ii)+439/216*Y_1-8*Y_2+3680/513*Y_3-845/4104*Y_4               ,m,k,s,f0,w);
    Y_6      = h*fvbr(t+1/2*h, y(:,ii)+8/27*Y_1-2*Y_2+3544/2565*Y_3+1859/4104*Y_4-11/40*Y_5  ,m,k,s,f0,w);
        
    % 4th order approximation:
    sol4(:,ii+1) = y(:,ii) + 25/216*Y_1+1408/2565*Y_3+2197/4101*Y_4-1/5*Y_5;
    % 5th order approximation:
    sol5(:,ii+1) = y(:,ii) + 16/135*Y_1+6656/12.825*Y_3+28.561/56.430*Y_4-9/50*Y_5+2/55*Y_6;
    % Error in solution (for y(:,ii+1) between 4th & 5th order RK.)
    R = abs(max(sol5(:,ii+1)-sol4(:,ii+1)))/h;   
    delta = 0.84*(epsilon/R)^(1/4);
    
    if R<=epsilon
        t = t+h;
        y(:,ii+1) = sol4(:,ii+1);
        ii=ii+1;
%         h = delta*h;
    else 
        h = delta*h;
        Y_1      = h*fvbr(t,y(:,ii)                                                              ,m,k,s,f0,w);
        Y_2      = h*fvbr(t+1/4*h, y(:,ii)+1/4*Y_1                                               ,m,k,s,f0,w);
        Y_3      = h*fvbr(t+3/8*h, y(:,ii)+3/32*Y_1+9/32*Y_2            ,m,k,s,f0,w);
        Y_4      = h*fvbr(t+12/13*h, y(:,ii)+1932/2197*Y_1-7200/2197*Y_2+7296/2197*Y_3           ,m,k,s,f0,w);
        Y_5      = h*fvbr(t+h, y(:,ii)+439/216*Y_1-8*Y_2+3680/513*Y_3-845/4104*Y_4               ,m,k,s,f0,w);
%       Y_6      = h*fvbr(t+1/2*h, y(:,ii)+8/27*Y_1-2*Y_2+3544/2565*Y_3+1859/4104*Y_4-11/40*Y_5  ,m,k,s,f0,w);
        y(:,ii+1) = y(:,ii) + 25/216*Y_1+1408/2565*Y_3+2197/4101*Y_4-1/5*Y_5;
        ii=ii+1;
    end
end


figure
subplot(1,2,1)
plot(t,y(1,:),'red')
xlabel('Time')
subplot(1,2,2)
plot(y(1,:),y(2,:))
title('Phase portrait')
y(1,end)

end