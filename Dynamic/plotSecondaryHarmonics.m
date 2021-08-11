% Plotting harmonics of secondary field

close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 1;
beta = 0.25;

M = 80;
rhos = linspace(1e-3, 3, M);
phi = linspace(-pi, pi, 100);
kz = linspace(-4, 4, 10);
omega = linspace(-2, 2, 10);
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

for n=0
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
            - 1*EzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        EphiFourier = 1*EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - 1*EphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        ErhoFourier = 1*ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - 1*ErhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        eta0HzFourier = 1*eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - 1*eta0HzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        eta0HphiFourier = 1*eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - 1*eta0HphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        eta0HrhoFourier = 1*eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - 1*eta0HrhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
    
        EzFFT = fftshift(fft2(fftshift(EzFourier))) * dkz * domega;
        EphiFFT = fftshift(fft2(fftshift(EphiFourier))) * dkz * domega;
        ErhoFFT = fftshift(fft2(fftshift(ErhoFourier))) * dkz * domega;
        eta0HzFFT = fftshift(fft2(fftshift(eta0HzFourier))) * dkz * domega;
        eta0HphiFFT = fftshift(fft2(fftshift(eta0HphiFourier))) * dkz * domega;
        eta0HrhoFFT = fftshift(fft2(fftshift(eta0HrhoFourier))) * dkz * domega;
        
%         EzTotal(i) = EzTotal(i) + EzFFT(5, 5);
%         EphiTotal(i) = EphiTotal(i) + EphiFFT(5, 5);
%         ErhoTotal(i) = ErhoTotal(i) + ErhoFFT(5, 5);
%         eta0HzTotal(i) = eta0HzTotal(i) + eta0HzFFT(5, 5);
%         eta0HphiTotal(i) = eta0HphiTotal(i) + eta0HphiFFT(5, 5);
%         eta0HrhoTotal(i) = eta0HrhoTotal(i) + eta0HrhoFFT(5, 5);
        
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