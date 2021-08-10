% Checking Bessel function sum convergence rate.

% Clearing workspace
close all;
clear;
clc;

beta = 0.05;
gamma = (1 - beta^2) ^ (-0.5);
rho = 0.5;
n = 0;

M = 1e2;
kz = 0.2;
omega = 0.2;
[Kz, W] = meshgrid(kz, omega);

omegaNorm = W / (gamma*beta);
hypot = sqrt(Kz.^2 + omegaNorm.^2);   

% Calculating sum w. maximal number of terms
nuMax = 80;
besselSumMax = 0;

tmp = zeros(1, 161);
i = 1;

tic
for nu=(-nuMax:nuMax)
    besselSumMax = besselSumMax + besselj(nu, (W/beta) .* rho) .* besseli(-n-nu, hypot.*rho);
    disp(besselSumMax);
    tmp(i) = besselj(nu, (W/beta) .* rho) .* besseli(-n-nu, hypot.*rho);
    i = i + 1;
end
toc

figure; hold on;
plot(log10(abs(tmp)), 'LineWidth', 1);

besselSumMaxNorm = norm(besselSumMax, 'fro');
disp(besselSumMaxNorm);

besselSumTail = besselSumMax;
nuVec = 0:80;
errVec = zeros(1, numel(nuVec));

for i=1:numel(nuVec)
    nu_i = nuVec(i);
    disp(nu_i);
    
    if nu_i == 0
        besselSumTail = besselSumTail - besselj(nu_i, (W/beta) .* rho) .* besseli(-n-nu_i, hypot.*rho);
    else
        for nu=[-nu_i, nu_i]
            besselSumTail = besselSumTail - besselj(nu, (W/beta) .* rho) .* besseli(-n-nu, hypot.*rho);
        end
    end
    errVec(i) = norm(besselSumTail, 'fro') / besselSumMaxNorm;
end

figure; hold on;
plot(nuVec, log10(errVec), 'LineWidth', 1);

    
