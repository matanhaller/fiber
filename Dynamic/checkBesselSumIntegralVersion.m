% Checking integral version of Bessel function sum convergence rate.

% Clearing workspace
close all;
clear;
clc;

beta = 0.99;
gamma = (1 - beta^2) ^ (-0.5);
rho = 1;
y0 = 1.4;

M = 1e2;
kz = linspace(-20, 20, 2e2);
omega = linspace(-20, 20, 2e2);
[Kz, W] = meshgrid(kz, omega);

omegaNorm = W / (gamma*beta);
hypot = sqrt(Kz.^2 + omegaNorm.^2);   

nuMax = 1e3;
nuVec = 1:nuMax;

figure;

tic
for n=0:4:20
    % Calculating sum w. maximal number of terms
    besselSumMax = exp(-hypot*abs(y0)).*besselSum(rho, n, Kz, W, beta, nuMax, 0);
    besselSumMaxNorm = norm(besselSumMax, 'fro');

    errVec = zeros(1, numel(nuVec));
    
    for i=1:numel(nuVec)
        disp(i);
        nu = nuVec(i);
        besselSumCurr = exp(-hypot*abs(y0)).*besselSum(rho, n, Kz, W, beta, nu, 0);
        errVec(i) = norm(besselSumCurr - besselSumMax, 'fro') / besselSumMaxNorm;
        if errVec(i) < 1e-6
            break;
        end
    end

    semilogy(nuVec, errVec, 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n)); hold on;
    xlabel('$N$', 'FontSize', 14);
    ylabel('Error', 'FontSize', 14);
end
toc

legend('FontSize', 14);
