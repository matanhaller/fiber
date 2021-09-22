% Plotting force exerted on the point charge by the fiber

close all;
clear;
clc;

M = 10;

x0 = 0; y0 = logspace(log10(1.02), 1, M); z0 = 0;
beta = 0.2;

% kz = logspace(-3, 1, 200); kz = [-flip(kz), kz];
% omega = logspace(log10(1.0001e-3), log10(10.0001), 200); omega = [-flip(omega), omega];
kz = linspace(-10, 10, 40);
omega = linspace(-10.0001, 10.0001, 40);
dkz = kz(2) - kz(1);
domega= omega(2) - omega(1);
[K, W] = meshgrid(kz, omega);
N = -10:10;
nuMax = 40;

%% Plotting force as function of distance from cylinder
figure; hold on;

tic
for epsilon=1.2
    
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

            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, rho, z0, beta, epsilon, nuMax);
            EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
            EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            
            EzTotal = trapz(omega, trapz(kz, EzFourier, 2));
            EphiTotal = trapz(omega, trapz(kz, EphiFourier, 2));
            ErhoTotal = trapz(omega, trapz(kz, ErhoFourier, 2));
            eta0HzTotal = trapz(omega, trapz(kz, eta0HzFourier, 2));
            eta0HphiTotal = trapz(omega, trapz(kz, eta0HphiFourier, 2));
            eta0HrhoTotal = trapz(omega, trapz(kz, eta0HrhoFourier, 2));

            EzInv(i) = EzInv(i) + (-1j)^n*EzTotal;
            EphiInv(i) = EphiInv(i) + (-1j)^n*EphiTotal;
            ErhoInv(i) = ErhoInv(i) + (-1j)^n*ErhoTotal;
            eta0HzInv(i) = eta0HzInv(i) + (-1j)^n*eta0HzTotal;
            eta0HphiInv(i) = eta0HphiInv(i) + (-1j)^n*eta0HphiTotal;
            eta0HrhoInv(i) = eta0HrhoInv(i) + (-1j)^n*eta0HrhoTotal;
        end
    end

    F = sqrt(real(EzInv + beta * eta0HrhoInv).^2 + real(EphiInv).^2 + real(ErhoInv - beta * eta0HzInv).^2);
    F0 = 0.25 * (epsilon - 1) / (epsilon + 1) * 1 ./ ((y0 - 1) .^ 2);

    plot(y0-1, F./F0, 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%.1f$', epsilon));
    xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$|\bar{F}|$', 'FontSize', 14, 'Interpreter', 'latex');
end
toc

legend('FontSize', 14, 'Interpreter', 'latex');

%% Plotting force as function of velocity
rho = 1.4;

figure(1); hold on;
figure(2); hold on;

tic
for epsilon=1.2
    betaC = 1 / sqrt(epsilon);
    gbC = betaC / sqrt(1 - betaC^2);
    gb = (0.1:0.1:4) * gbC; 
    betas = sqrt(gb.^2 ./ (gb.^2 + 1));
    Fr = zeros(1, numel(betas));
    Fa = zeros(1, numel(betas));
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

            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, rho, z0, beta, epsilon, nuMax);
            EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
            EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);            
            ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);

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

        Fr(i) = -real(ErhoInv - beta * eta0HzInv);
        disp(real(ErhoInv));
        disp(real(eta0HzInv));
        Fa(i) = -real(EphiInv);
    end
    
    figure(1);
    plot(gb./gbC, Fr, 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%.1f$', epsilon));
    figure(2);
    plot(gb./gbC, Fa, 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%.1f$', epsilon));
end
toc

for i=1:2
    figure(i);
    xlabel('$\gamma \beta / (\gamma \beta)_\mathrm{C}$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$|\bar{F}|$', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end
