% Checking relative error in field Fourier components as function of the
% truncation limit.

close all;
clear;
clc;

epsilon = 4;
beta = 0.8;
x0 = 0; y0 = 1.4; z0 = 0;
rho = 0.75;
N = 0:20;
kz = linspace(-10, 10, 200);
omega = linspace(-10.0001, 10.0001, 200);
[K, W] = meshgrid(kz, omega);

nuMaxIdeal = 200;
nuVec = 0:45;

tic;

for i=1:6
    figure(i); hold on;
end

for n=10
    disp(n);
    [AnkIdeal, BnkIdeal, eta0CnkIdeal, eta0DnkIdeal] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMaxIdeal);

    EzIdeal = EzSecondaryFourier(AnkIdeal, BnkIdeal, rho, n, K, W, epsilon);
    EphiIdeal = EphiSecondaryFourier(AnkIdeal, BnkIdeal, eta0CnkIdeal, eta0DnkIdeal, rho, n, K, W, epsilon);
    ErhoIdeal = ErhoSecondaryFourier(AnkIdeal, BnkIdeal, eta0CnkIdeal, eta0DnkIdeal, rho, n, K, W, epsilon);
    eta0HzIdeal = eta0HzSecondaryFourier(eta0CnkIdeal, eta0DnkIdeal, rho, n, K, W, epsilon);
    eta0HphiIdeal = eta0HphiSecondaryFourier(AnkIdeal, BnkIdeal, eta0CnkIdeal, eta0DnkIdeal, rho, n, K, W, epsilon);
    eta0HrhoIdeal = eta0HrhoSecondaryFourier(AnkIdeal, BnkIdeal, eta0CnkIdeal, eta0DnkIdeal, rho, n, K, W, epsilon);

    EzIdealNorm = norm(EzIdeal, 'fro');
    EphiIdealNorm = norm(EphiIdeal, 'fro');
    ErhoIdealNorm = norm(ErhoIdeal, 'fro');
    eta0HzIdealNorm = norm(eta0HzIdeal, 'fro');
    eta0HphiIdealNorm = norm(eta0HphiIdeal, 'fro');
    eta0HrhoIdealNorm = norm(eta0HrhoIdeal, 'fro');

    EzRMSE = zeros(1, numel(nuVec));
    EphiRMSE = zeros(1, numel(nuVec));
    ErhoRMSE = zeros(1, numel(nuVec));
    eta0HzRMSE = zeros(1, numel(nuVec));
    eta0HphiRMSE = zeros(1, numel(nuVec));
    eta0HrhoRMSE = zeros(1, numel(nuVec));

    for i=1:numel(nuVec)
        nu = nuVec(i);
        disp(nu);

        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nu);

        Ez = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
        Ephi = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        Erho = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0Hz = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0Hphi = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0Hrho = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);

        EzRMSE(i) = norm(Ez - EzIdeal, 'fro') / EzIdealNorm;
        EphiRMSE(i) = norm(Ephi - EphiIdeal, 'fro') / EphiIdealNorm;
        ErhoRMSE(i) = norm(Erho - ErhoIdeal, 'fro') / ErhoIdealNorm;
        eta0HzRMSE(i) = norm(eta0Hz - eta0HzIdeal, 'fro') / eta0HzIdealNorm;
        eta0HphiRMSE(i) = norm(eta0Hphi - eta0HphiIdeal, 'fro') / eta0HphiIdealNorm;
        eta0HrhoRMSE(i) = norm(eta0Hrho - eta0HrhoIdeal, 'fro') / eta0HrhoIdealNorm;
    end
    
    figure(1);
    plot(nuVec, log10(EzRMSE), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(2);
    plot(nuVec, log10(EphiRMSE), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(3);
    plot(nuVec, log10(ErhoRMSE), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(4);
    plot(nuVec, log10(eta0HzRMSE), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(5);
    plot(nuVec, log10(eta0HphiRMSE), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(6);
    plot(nuVec, log10(eta0HrhoRMSE), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));

end

figure(1);
title('Longitudinal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(2);
title('Azimuthal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(3);
title('Radial Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(4);
title('Longitudinal Magnetic Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(5);
title('Azimuthal Magnetic Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(6);
title('Radial Magnetic Field', 'FontSize', 14, 'Interpreter', 'latex');

for i=1:6
    figure(i);
    xlabel('$\nu$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Relative Error', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end

toc;