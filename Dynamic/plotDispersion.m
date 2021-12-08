% Plotting dispersion curve of the dielectric cylinder.

close all;
clear;
clc;

epsilon = 2;
N = -10:10;

M = 2.5e3;
maxVal = 10;
maxModes = 15;
omega = linspace(1e-3, maxVal+1e-4, M);

tic;

for n=3
    kz = zeros(maxModes, numel(omega));
    
    for i=1:numel(omega)
        w = omega(i);
        kRange = linspace(w, sqrt(epsilon)*w-1e-4, 1e4); kRange(1) = [];
        d = DeltaCylinder(n, kRange, w, epsilon);
        sp = spmak(augknt(kRange,2), real(d));
        k0 = fnzeros(sp);
        kz(1:size(k0,2),i) = k0(1,:).';
    end

    kz = sort(kz, 1, 'descend');
    
    % Fixing 1st eigenmode for n=+-1
    if abs(n) == 1
        kz(1,kz(1,:) == 0) = omega(kz(1,:) == 0) + 1e-4;
    end

    figure; hold on;
    for i=1:size(kz,1)
        kzi = kz(i,:);
        plot(omega(kzi > 0), kzi(kzi > 0), 'LineWidth', 2);
    end

    plot(omega, omega, 'k--', 'LineWidth', 1);
    plot(omega, omega*sqrt(epsilon), 'k--', 'LineWidth', 1);

    xlim([0, maxVal]);
    ylim([0, maxVal+10]);

    xlabel('$\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$\bar{k}_z$', 'FontSize', 14, 'Interpreter', 'latex');
    
    kz = kz(1:maxModes,:);
%     save(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, n), 'kz');
end

toc;