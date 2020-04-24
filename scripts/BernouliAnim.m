h = figure;


BB = load('Bernoulli_MC_2_bigeps.mat');
Montecarlos = BB.Montecarlos;
yfref = BB.yfref;
time2 = BB.time2;
ndof = BB.ndof;
dt = pi/101;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Bernoulli_perturb.gif';
y = Montecarlos(:,:, 3000);
for i = 1:101
    clf(h,'reset')
    befores = y(i, 1:2:ndof)+ 270*y(i, 2:2:ndof);
    beforesX = 1:1:ndof/2;
    beforesX = beforesX - 0.1;
    afterX = 1:1:ndof/2;
    afterX = afterX + 0.1;
    
    afters = y(i, 1:2:ndof) - 270*y(i, 2:2:ndof);
    
    hold on;
    plot(y(i, 1:2:ndof), 'o', 'Linewidth', 5, 'color', 'blue');
    plot(beforesX, befores, '-o', 'Linewidth', 5, 'color', 'blue');
    plot(afterX, afters, '-o', 'Linewidth', 5, 'color', 'blue');
    %plot(yfref(i, 1:2:ndof), '-o', 'Linewidth', 5, 'color', 'Black'); 
    xlim([1, 4]);
    ylim([-0.01, 0.35]);
    ylabel('w(x)');
    xlabel('x');
    hold off;
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',3*dt); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',3*dt); 
      end 
end