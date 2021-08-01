close all;
clear;
clc;

N = 20;
K = 100;

r = linspace(1e-3, 5, 200);
theta = linspace(0, 2*pi, 200);
[R, T] = meshgrid(r, theta);
X = R.*cos(T); Y = R.*sin(T);
z=0;

etas = flip(linspace(1.2,3.8,20));
epsilonR = 12;

k = figure;
set(gcf, 'color', 'white');
set(gcf, 'Renderer', 'Painter');
filename_1 = 'changingEtas.gif';

for i=1:numel(etas)
    eta = etas(i);
    phiP = primaryPotentialOfPointCharge(eta, X, Y, z);
    phiS = secondaryPotential(epsilonR, eta, X, Y, z, N, K, 1e-3);
    phiTotal = phiP + phiS;
    plot(0,eta,'r.', 'MarkerSize',20);
    hold on
    surf(X, Y, phiTotal, 'EdgeColor', 'none');
    
    plot(cos(theta), sin(theta), '--w', 'LineWidth', 2);
    
    view(2);
    colorbar;
    caxis([-0.75, -0.1]);
    hold off
    xlim([-3, 3]);
    ylim([-2, 4]);
   
    title(sprintf('$\\eta=%.3f$', eta), 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    drawnow 
      % Capture the plot as an image 
      frame = getframe(k); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename_1,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename_1,'gif','WriteMode','append'); 
      end 
end