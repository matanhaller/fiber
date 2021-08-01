% Plotting a 2D plot of the intensity in different cross-of the Cylinder

close all;
clear;
clc;

epsilonR = [1, 2, 4, 12];
eta = 1.8;
N = 20;
K = 100;

% x = linspace(-3, 3, 40); y = linspace(-3, 3, 40); z = 0;
r = linspace(1e-3, 1, 40);
theta = linspace(0, 2*pi, 100);
z = 0;
[R, T] = meshgrid(r, theta);
X = R.*cos(T); Y = R.*sin(T);
% [X, Y] = meshgrid(x, y);

figure;
set(gcf, 'Renderer', 'Painter');

for i=1:4
    er = epsilonR(i);
    int = intensity(er, eta, X, Y, z, N, K, 1e-3);
    
    subplot(2, 2, i); hold on;
    surf(X, Y, log10(int), 'EdgeColor', 'none');
    view(2);
    colorbar;
    caxis([-2.5, 0.5]);
    title(sprintf('$\\varepsilon=%d$', er), 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
%              saveas(gcf, 'C:\Users\yansp\Dropbox\Matan and Yaniv\Static\Graphs\Intensity.jpg');
end


   
