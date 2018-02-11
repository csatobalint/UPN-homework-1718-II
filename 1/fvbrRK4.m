function fvbrRK4(T,N,m,k,s,f0,w,c1,c2)

h=T/N; % stepsize
t=linspace(0,T,N+1); % grid
y=zeros(2,N+1); % numerical solution

% RK4:
y(:,1)=[c1 c2]';
for j=1:N
    Y_1      = y(:,j);
    Y_2      = y(:,j)   + h / 2 * fvbr(t(j)         , Y_1,m,k,s,f0,w);
    Y_3      = y(:,j)   + h / 2 * fvbr(t(j) + h / 2 , Y_2,m,k,s,f0,w);
    Y_4      = y(:,j)   + h     * fvbr(t(j) + h / 2 , Y_3,m,k,s,f0,w);
    y(:,j+1) = y(:,j) + h / 6 *(fvbr(t(j),Y_1,m,k,s,f0,w)+ 2*fvbr(t(j) + h / 2,Y_2,m,k,s,f0,w)...
            + 2 * fvbr(t(j) + h / 2,Y_3,m,k,s,f0,w) + fvbr(t(j+1),Y_4,m,k,s,f0,w));

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