%% Unsteady flows in pipe networks homework. Applying Lax-Wendroff scheme (and other methods) for gas outflow from a tank.
close all
clc

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for n = 1:1:5
    % Draw plot for y = x.^n
    subplot(2,2,1)
    plot(x, p_write(n,:))
    title('Pressure')
    subplot(2,2,2)
    plot(x, v_write(n,:))
    title('Velocity')
    subplot(2,2,3)
    %plot(x, p_write(n,:))
    title('Density')
    subplot(2,2,4)
    %plot(x, v_write(n,:))
    title('Temperature')
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
  end