% Plotting a 2D plot of the total potential in different cross-sections

close all;
clear;
clc;

epsilonR = [1, 2, 4, 12];
eta = 1.8;
N = 20;
K = 100;

x = linspace(-3, 3, 20);
y = linspace(-3, 3, 20);
z = 0;
[X, Y] = meshgrid(x, y);

theta = linspace(0, 2*pi, 100);

phiP = primaryPotentialOfPointCharge(eta, X, Y, z);

% Plotting for different values of epsilonR
figure(1);
tiledlayout(2, 2);
set(gcf, 'Renderer', 'Painter');

for i=1:4
    er = epsilonR(i);
    phiS = secondaryPotential(er, eta, X, Y, z, N, K, 1e-3);
    phiTotal = phiP + phiS;
    
    nexttile; hold on;
    surf(X, Y, phiTotal, 'EdgeColor', 'none');
    plot(cos(theta), sin(theta), 'w', 'LineWidth', 1);
    view(2);
    caxis([-0.75, -0.19]);
    xlim([-3, 3]);
    ylim([-3, 3]);
    title(sprintf('$\\varepsilon=%d$', er), 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
end

cb = colorbar;
cb.Layout.Tile = 'east';

% Plotting for different values of eta
x = linspace(-2, 2, 20);
y = linspace(-2, 6, 20);
z = 0;
[X, Y] = meshgrid(x, y);

etas = [4, 2, 1.5, 1.2];
epsilonR = 12;

figure(2);
tiledlayout(1, 4);
set(gcf, 'Renderer', 'Painter');

for i=1:numel(etas)
    eta = etas(i);
    phiP = primaryPotentialOfPointCharge(eta, X, Y, z);
    phiS = secondaryPotential(epsilonR, eta, X, Y, z, N, K, 1e-3);
    phiTotal = phiP + phiS;
    
    nexttile; hold on;
    surf(X, Y, phiTotal, 'EdgeColor', 'none');
    plot(cos(theta), sin(theta), '--w', 'LineWidth', 2);
    view(2);
    caxis([-0.75, -0.1]);
    xlim([-2, 2]);
    ylim([-2, 6]);
    title(sprintf('$\\eta=%.1f$', eta), 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
end

cb = colorbar;
cb.Layout.Tile = 'east';
