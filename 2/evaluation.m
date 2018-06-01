%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Run first the main.m                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting and saving results
h = figure('visible','off');
x0=50;
y0=50;
width=800;
height=600;
desiredFps=30;
desiredVideoLength=8;
set(h,'units','points','position',[x0,y0,width,height])
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'p-v-rho-T.gif';
filename2 = ['plots/p-v-rho-T_dx_' num2str(dx) ''];
vobj=VideoWriter(filename2, 'Motion JPEG AVI');
vobj.FrameRate=desiredFps;
vobj.Quality=100;
open(vobj);
saveRate=round(last_timestep/(desiredVideoLength*desiredFps));
fontSizeLabels=14;
for n = 1:saveRate:last_timestep
    suptitle(['Time: ' num2str(t(n)) ' [s] Grid size: ' num2str(dx) '[m]'])
    subplot(2,2,1)
        area(x, p_write(n,:))
        title(' ')
        xlim([xL xR])
        xlabel('$x~[m]$','interpreter','latex')
        ylim([p_ref*0.95 p_res*1.1])
        ylabel('$p~[Pa]$','interpreter','latex')
        set(gca,'FontSize',fontSizeLabels)
    subplot(2,2,2)
        area(x, v_write(n,:))
        title(' ')
        xlim([xL xR])
        xlabel('$x~[m]$','interpreter','latex')
        ylim([0 1.5*v_write(last_timestep,end)])
        %ylim([-10 10]) %if wall are BCs
        ylabel('$v~[\frac{m}{s}]$','interpreter','latex')
        set(gca,'FontSize',fontSizeLabels)
    subplot(2,2,3)
        area(x, rho_write(n,:))
        %title('Density')
        xlim([xL xR])
        xlabel('$x~[m]$','interpreter','latex')
        ylim([1 1.5*rho_write(last_timestep,end)])
        ylabel('$\rho~[\frac{kg}{m^3}]$','interpreter','latex')
        set(gca,'FontSize',fontSizeLabels)    
    subplot(2,2,4)
        area(x, T_write(n,:))
        %title('Temperature')
        xlim([xL xR])
        xlabel('$x~[m]$','interpreter','latex')
        ylim([250 350])
        ylabel('$T~[K]$','interpreter','latex')
        set(gca,'FontSize',fontSizeLabels)
%     drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif','Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
      
     writeVideo(vobj, frame);
     cla(gca)
end
 close(vobj)
toc

%% Determining the final (stationary) Mach number:
% t_quasistat=floor(0.3*size(find(t),2));
% M_st=v_write(t_quasistat,:)./sqrt(gamma*R*T_write(t_quasistat,:));
M_st=v_write(size(find(t),2),:)./sqrt(gamma*R*T_write(size(find(t),2),:));
M_st_av=mean(M_st);
close all
M_cso=KiaramlasFanno(gamma, T(1)-273, p_res/10^5, p_ref/10^5, L, d*1000, lambda, R,dx)
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(x,M_cso,'DisplayName','Numerical solution','MarkerFaceColor',[1 0 0],...
    'MarkerSize',6,...
    'Marker','o',...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',[1 0 0]);
plot(x,M_st,'DisplayName','Fanno Flow','MarkerSize',8,'Marker','x',...
    'LineWidth',1,...
    'Color',[0 0 1]);

xlabel('$x~[m]$','Interpreter','latex');
ylabel('$Ma~[-]$','Interpreter','latex');
titlename=['Grid size: ' num2str(dx) ' [m]'];
title(titlename);
xlim([xL xR]);
ylim([M_st(1)*0.9 M_st(end)*1.1]);
box(axes1,'on');
grid(axes1,'on');
set(axes1,'FontSize',14);
legend1 = legend(axes1,'show');
filename = ['plots/ma_dx_' num2str(dx) '.png'];
saveas(gcf,filename)

%% Plotting and saving results
% WALL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
h = figure('visible','off');
x0=50;
y0=50;
width=800;
height=600;
desiredFps=30;
desiredVideoLength=8;
set(h,'units','points','position',[x0,y0,width,height])
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'p-v-rho-T.gif';
filename2 = ['plots/p-v-rho-T_dx_' num2str(dx) ''];
filename2 = ['plots/wall_2_dx' num2str(dx) ''];
vobj=VideoWriter(filename2, 'Motion JPEG AVI');
vobj.FrameRate=desiredFps;
vobj.Quality=100;
open(vobj);
saveRate=round(last_timestep/(desiredVideoLength*desiredFps));
fontSizeLabels=14;
for n = 1:saveRate:last_timestep
    suptitle(['Time: ' num2str(t(n)) ' [s] Grid size: ' num2str(dx) '[m]'])
    subplot(2,1,1)
        area(x, p_write(n,:))
        title(' ')
        xlim([xL xR])
        xlabel('$x~[m]$','interpreter','latex')
        ylim([0.9*10^5 1.3*10^5])
        ylabel('$p~[Pa]$','interpreter','latex')
        set(gca,'FontSize',fontSizeLabels)
    subplot(2,1,2)
        area(x, v_write(n,:))
        title(' ')
        xlim([xL xR])
        xlabel('$x~[m]$','interpreter','latex')
        ylim([-30 30])
        ylabel('$v~[\frac{m}{s}]$','interpreter','latex')
        set(gca,'FontSize',fontSizeLabels)
%     drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif','Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
      
     writeVideo(vobj, frame);
     cla(gca)
end
 close(vobj)
toc

