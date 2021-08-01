% Plotting harmonics of secondary field

close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 1;
beta = 0.25;

M = 40;
rhos = linspace(1e-3, 0.95, M);
phi = linspace(-pi, pi, 100);
kz = linspace(-6.1, 6.1, 40);
omega = linspace(-4.1, 4.1, 40);
[K, W] = meshgrid(kz, omega);

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);

[R, P] = meshgrid(rhos, phi);
X = R.*cos(P); Y = R.*sin(P);

EzInv = zeros(size(R));
EphiInv = zeros(size(R));
ErhoInv = zeros(size(R));
eta0HzInv = zeros(size(R));
eta0HphiInv = zeros(size(R));
eta0HrhoInv = zeros(size(R));

for n=-2:2
    EzTotal = zeros(1, numel(rhos));
    EphiTotal = zeros(1, numel(rhos));
    ErhoTotal = zeros(1, numel(rhos));
    eta0HzTotal = zeros(1, numel(rhos));
    eta0HphiTotal = zeros(1, numel(rhos));
    eta0HrhoTotal = zeros(1, numel(rhos));

    for i=1:numel(rhos)
        rho = rhos(i);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon);
        EzFourier = 1*EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            + 1*EzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(rho - 1);
        EphiFourier = 1*EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            + 1*EphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(rho - 1);
        ErhoFourier = 1*ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            + 1*ErhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(rho - 1);
        eta0HzFourier = 1*eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            + 1*eta0HzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(rho - 1);
        eta0HphiFourier = 1*eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            + 1*eta0HphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(rho - 1);
        eta0HrhoFourier = 1*eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            + 1*eta0HrhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(rho - 1);

        EzTotal(i) = EzTotal(i) + sum(sum(EzFourier)) * dkz * domega;
        EphiTotal(i) = EphiTotal(i) + sum(sum(EphiFourier)) * dkz * domega;
        ErhoTotal(i) = ErhoTotal(i) + sum(sum(ErhoFourier)) * dkz * domega;
        eta0HzTotal(i) = eta0HzTotal(i) + sum(sum(eta0HzFourier)) * dkz * domega;
        eta0HphiTotal(i) = eta0HphiTotal(i) + sum(sum(eta0HphiFourier)) * dkz * domega;
        eta0HrhoTotal(i) = eta0HrhoTotal(i) + sum(sum(eta0HrhoFourier)) * dkz * domega;

        disp(rho);
    end
    EzInv = EzInv + exp(1j*n*phi)'*EzTotal;
    EphiInv = EphiInv + exp(1j*n*phi)'*EphiTotal;
    ErhoInv = ErhoInv + exp(1j*n*phi)'*ErhoTotal;
    eta0HzInv = eta0HzInv + exp(1j*n*phi)'*eta0HzTotal;
    eta0HphiInv = eta0HphiInv + exp(1j*n*phi)'*eta0HphiTotal;
    eta0HrhoInv = eta0HrhoInv + exp(1j*n*phi)'*eta0HrhoTotal;
end

I = abs(EzInv).^2 + abs(EphiInv).^2 + abs(ErhoInv).^2 ...
  + abs(eta0HzInv).^2 + abs(eta0HphiInv).^2 + abs(eta0HrhoInv).^2;

figure; hold on;
surf(X, Y, I, 'EdgeColor', 'none');
view(2);
xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
title(sprintf('Longitudinal Field Intensity: $N=%d$, $\\varepsilon=%.2f$', 0, epsilon), 'FontSize', 14, 'Interpreter', 'latex');
colorbar;