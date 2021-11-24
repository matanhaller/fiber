% Checking Bessel function sum convergence rate.

% Clearing workspace
close all;
clear;
clc;

beta = 0.2;
gamma = (1 - beta^2) ^ (-0.5);
rho = 1;
y0 = 1.4;

M = 1e2;
kz = linspace(-10, 10, 1e2);
omega = linspace(-10, 10, 1e2);
[Kz, W] = meshgrid(kz, omega);

omegaNorm = W / (gamma*beta);
hypot = sqrt(Kz.^2 + omegaNorm.^2);   

figure;
for n=0:4:40
    % Calculating sum w. maximal number of terms
    nuMax = 200;
    besselSumMax = 0;

    tic
    for nu=(-nuMax:nuMax)
        besselSumMax = besselSumMax + exp(-hypot*abs(y0)).*(besselj(nu, (W/beta) .* rho) .* besseli(-n-nu, hypot.*rho));
    end
    toc

    besselSumMaxNorm = norm(besselSumMax, 'fro');

    besselSumTail = besselSumMax;
    nuVec = 0:200;
    errVec = zeros(1, numel(nuVec));

    for i=1:numel(nuVec)
        nu_i = nuVec(i);
        disp(nu_i);

        if nu_i == 0
            besselSumTail = besselSumTail - exp(-hypot*abs(y0)).*(besselj(nu_i, (W/beta) .* rho) .* besseli(-n-nu_i, hypot.*rho));
        else
            for nu=[-nu_i, nu_i]
                besselSumTail = besselSumTail - exp(-hypot*abs(y0)).*(besselj(nu, (W/beta) .* rho) .* besseli(-n-nu, hypot.*rho));
            end
        end
        errVec(i) = norm(besselSumTail, 'fro') / besselSumMaxNorm;
    end

    semilogy(nuVec, errVec, 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n)); hold on;
    xlabel('$\nu$', 'FontSize', 14);
    ylabel('Error', 'FontSize', 14);
end

legend('FontSize', 14);
