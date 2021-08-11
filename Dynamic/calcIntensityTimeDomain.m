% Calculating intesntiy in the time domain and saving workspace.

% Clearing workspace
close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 0.8;

M = 120;
rhos = linspace(1e-3, 3, M);
kz = linspace(-8, 8, 40);
omega = linspace(-8.0001, 8.0001, 40);
[K, W] = meshgrid(kz, omega);

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);

z = linspace(-pi/dkz, pi/dkz, numel(kz)+1); z(end) = [];
t = linspace(-pi/domega, pi/domega, numel(omega)+1); t(end) = [];

tic
for n=-2:2
    disp(n);
    
    EzTotal = zeros(numel(omega), numel(kz), numel(rhos));
    EphiTotal = zeros(numel(omega), numel(kz), numel(rhos));
    ErhoTotal = zeros(numel(omega), numel(kz), numel(rhos));
    eta0HzTotal = zeros(numel(omega), numel(kz), numel(rhos));
    eta0HphiTotal = zeros(numel(omega), numel(kz),numel(rhos));
    eta0HrhoTotal = zeros(numel(omega), numel(kz), numel(rhos));

    for i=1:numel(rhos)
        rho = rhos(i);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon);
        
        EzFourier = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - EzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - EphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - ErhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        eta0HzFourier = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - eta0HzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - eta0HphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, x0, y0, z0, beta, epsilon) ...
            - eta0HrhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta) * heaviside(1 - rho);
        
        EzIFFT = fftshift(fft2(fftshift(EzFourier))) * dkz * domega;
        EphiIFFT = fftshift(fft2(fftshift(EphiFourier))) * dkz * domega;
        ErhoIFFT = fftshift(fft2(fftshift(ErhoFourier))) * dkz * domega;
        eta0HzIFFT = fftshift(fft2(fftshift(eta0HzFourier))) * dkz * domega;
        eta0HphiIFFT = fftshift(fft2(fftshift(eta0HphiFourier))) * dkz * domega;
        eta0HrhoIFFT = fftshift(fft2(fftshift(eta0HrhoFourier))) * dkz * domega;
        
        EzTotal(:,:,i) = flip(EzIFFT, 1);
        EphiTotal(:,:,i) = flip(EphiIFFT, 1);
        ErhoTotal(:,:,i) = flip(ErhoIFFT, 1);
        eta0HzTotal(:,:,i) = flip(eta0HzIFFT, 1);
        eta0HphiTotal(:,:,i) = flip(eta0HphiIFFT, 1);
        eta0HrhoTotal(:,:,i) = flip(eta0HrhoIFFT, 1);
        
        disp(rho);
    end
    save(sprintf('Time Domain/Ez_n=%d.mat', n), 'EzTotal');
    save(sprintf('Time Domain/Ephi_n=%d.mat', n), 'EphiTotal');
    save(sprintf('Time Domain/Erho_n=%d.mat', n), 'ErhoTotal');
    save(sprintf('Time Domain/eta0Hz_n=%d.mat', n), 'eta0HzTotal');
    save(sprintf('Time Domain/eta0Hphi_n=%d.mat', n), 'eta0HphiTotal');
    save(sprintf('Time Domain/eta0Hrho_n=%d.mat', n), 'eta0HrhoTotal');
end
toc
