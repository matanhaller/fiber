% Plotting force exerted on the point charge by the fiber

close all;
clear;
clc;

M = 10;

x0 = 0; y0 = logspace(log10(1.02), 1, M); z0 = 0;
beta = 0.1;

kz = linspace(-8, 8, 200);
omega = linspace(-6, 6, 200);
[K, W] = meshgrid(kz, omega);
N = -10:10;

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);

%% Plotting force as function of distance from cylinder
figure; hold on;

tic
for epsilon=[1.2, 2, 4, 12]
    
    EzInv = zeros(1, M);
    EphiInv = zeros(1, M);
    ErhoInv = zeros(1, M);
    eta0HzInv = zeros(1, M);
    eta0HphiInv = zeros(1, M);
    eta0HrhoInv = zeros(1, M);
    
    for i=1:M
        rho = y0(i);
        disp(rho);

        for n=N
            disp(n);

            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, rho, z0, beta, epsilon);
            EzFourier = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            eta0HzFourier = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            
            EzTotal = sum(sum(EzFourier)) * dkz * domega;
            EphiTotal = sum(sum(EphiFourier)) * dkz * domega;
            ErhoTotal = sum(sum(ErhoFourier)) * dkz * domega;
            eta0HzTotal = sum(sum(eta0HzFourier)) * dkz * domega;
            eta0HphiTotal = sum(sum(eta0HphiFourier)) * dkz * domega;
            eta0HrhoTotal = sum(sum(eta0HrhoFourier)) * dkz * domega;

            EzInv(i) = EzInv(i) + (-1j)^n*EzTotal;
            EphiInv(i) = EphiInv(i) + (-1j)^n*EphiTotal;
            ErhoInv(i) = ErhoInv(i) + (-1j)^n*ErhoTotal;
            eta0HzInv(i) = eta0HzInv(i) + (-1j)^n*eta0HzTotal;
            eta0HphiInv(i) = eta0HphiInv(i) + (-1j)^n*eta0HphiTotal;
            eta0HrhoInv(i) = eta0HrhoInv(i) + (-1j)^n*eta0HrhoTotal;
        end
    end

    F = sqrt(abs(EzInv + eta0HrhoInv).^2 + abs(EphiInv).^2 + abs(ErhoInv - eta0HzInv).^2);
    F0 = 0.25 * (epsilon - 1) / (epsilon + 1) * 1 ./ ((y0 - 1) .^ 2);

    plot(y0-1, F./F0, 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%.1f$', epsilon));
    xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$|\bar{F}|$', 'FontSize', 14, 'Interpreter', 'latex');
end
toc

legend('FontSize', 14, 'Interpreter', 'latex');

%% Plotting force as function of velocity
rho = 1.5;
M = 4;
F = zeros(1, M);

figure; hold on;

tic
for epsilon=[1.2, 2, 4, 12]
    betaC = 1 / sqrt(epsilon);
    gbC = betaC / sqrt(1 - betaC^2);
    gb = [0.5, 1, 2, 4] * gbC; 
    betas = sqrt(gb.^2 ./ (gb.^2 + 1));
    for i=1:numel(betas)
        beta = betas(i);
        disp(beta);

        EzInv = 0;
        EphiInv = 0;
        ErhoInv = 0;
        eta0HzInv = 0;
        eta0HphiInv = 0;
        eta0HrhoInv = 0;

        for n=N
            disp(n);

            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, rho, z0, beta, epsilon);
            EzFourier = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);            
            ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            eta0HzFourier = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);
            eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, rho, z0, beta, epsilon);

            EzTotal = sum(sum(EzFourier)) * dkz * domega;
            EphiTotal = sum(sum(EphiFourier)) * dkz * domega;
            ErhoTotal = sum(sum(ErhoFourier)) * dkz * domega;
            eta0HzTotal = sum(sum(eta0HzFourier)) * dkz * domega;
            eta0HphiTotal = sum(sum(eta0HphiFourier)) * dkz * domega;
            eta0HrhoTotal = sum(sum(eta0HrhoFourier)) * dkz * domega;

            EzInv = EzInv + (-1j)^n*EzTotal;
            EphiInv = EphiInv + (-1j)^n*EphiTotal;
            ErhoInv = ErhoInv + (-1j)^n*ErhoTotal;
            eta0HzInv = eta0HzInv + (-1j)^n*eta0HzTotal;
            eta0HphiInv = eta0HphiInv + (-1j)^n*eta0HphiTotal;
            eta0HrhoInv = eta0HrhoInv + (-1j)^n*eta0HrhoTotal;
        end

        F(i) = sqrt(abs(EzInv + eta0HrhoInv).^2 + abs(EphiInv).^2 + abs(ErhoInv - eta0HzInv).^2);
    end
    plot(gb./gbC, F/F(1), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%.1f$', epsilon));
end
toc

xlabel('$\gamma \beta / (\gamma \beta)_\mathrm{C}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$|\bar{F}|$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

