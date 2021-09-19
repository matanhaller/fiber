% Plotting the spectrum of radiation due to propagating waves (in the far
% field regime).

close all;
clear;
clc;

N = -10:10;
epsilon = [2, 3, 4, 6, 8, 12];
x0 = 0; y0 = 1.5; z0 = 0;
beta = 0.49;

M = 1e3;
rho = 1.5;
kz = linspace(1e-4, 0.75, M);
omegas = linspace(0.01, 10, M);
W = zeros(numel(epsilon), numel(omegas));

%% Plotting spectrum for different permittivities
tic
figure;
for j=1:numel(epsilon)
    er = epsilon(j);
    disp(er);
    for n=N
        disp(n);
        for i=1:numel(omegas)
            omega = omegas(i);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kRange, omega, x0, y0, z0, beta, er);
            W(j, i) = W(j, i) + omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        end
    end
end

Wmax = max(max(W));
for j=1:numel(epsilon)
    er = epsilon(j);
    semilogy(omegas, W(j, :) / Wmax, 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%d$', er)); hold on;
end

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');
toc

%% Plotting spectrum for different velocities
epsilon = 4;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = [0.5, 0.9, 1.1, 2] * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
W = zeros(numel(betas), numel(omegas));

tic
figure; hold on;
for j=1:numel(betas)
    beta = betas(j);
    disp(beta);
    for n=N
        disp(n);
        for i=1:numel(omegas)
            omega = omegas(i);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kRange, omega, x0, y0, z0, beta, epsilon);
            W(j, i) = W(j, i) + omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        end
    end
end

for j=1:numel(betas)
    beta = betas(j);
    plot(omegas, W(j, :) / max(W(j, :)), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.2f$', beta));
end

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');
toc

%% Integrating over spectrum for different harmonics
epsilon = 12;
N = -10:10;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = [0.5, 1, 2, 4] * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
Energy = zeros(numel(betas), numel(N));

tic
figure;
for j=1:numel(betas)
    beta = betas(j);
    disp(beta);
    for k=1:numel(N)
        n = N(k);
        disp(n);
        W = zeros(1, numel(omegas));
        for i=1:numel(omegas)
            omega = omegas(i);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kRange, omega, x0, y0, z0, beta, epsilon);
            W(i) = omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        end
        Energy(j, k) = trapz(omegas, W);
    end
end

for j=1:numel(betas)
    beta = betas(j);
    semilogy(N, Energy(j, :) ./ max(Energy(j, :)), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.2f$', beta)); hold on;
end

toc

xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Energy', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');