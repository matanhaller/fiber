% Calculating intensity in the time domain and saving workspace.

% Clearing workspace
close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 0.8;
nuMax = 40;

M = 120;
rhos = linspace(1e-3, 5, M);
kz = linspace(-10, 10, 200);
% omega = [-2.822, 2.822];
omega = linspace(-10.0001, 10.0001, 200);
[K, W] = meshgrid(kz, omega);

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);
% domega = 1;

z = linspace(-pi/dkz, pi/dkz, numel(kz)+1); z(end) = [];
t = linspace(-pi/domega, pi/domega, numel(omega)+1); t(end) = [];

tic
for n=1
    disp(n);
    
    EzTotal = zeros(numel(omega), numel(kz), numel(rhos));
    EphiTotal = zeros(numel(omega), numel(kz), numel(rhos));
    ErhoTotal = zeros(numel(omega), numel(kz), numel(rhos));
    eta0HzTotal = zeros(numel(omega), numel(kz), numel(rhos));
    eta0HphiTotal = zeros(numel(omega), numel(kz),numel(rhos));
    eta0HrhoTotal = zeros(numel(omega), numel(kz), numel(rhos));
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax);
    
    for i=1:numel(rhos)
        rho = rhos(i);
        
        bs = besselSum(rho, n, K, W, beta, nuMax);
        bsDeriv = besselSumDeriv(rho, n, K, W, beta, nuMax);

        EzFourierP = EzPrimaryFourier(n, K, W, x0, y0, z0, beta, bs);
        EzFourierDerivP = EzPrimaryFourierDeriv(n, K, W, x0, y0, z0, beta, bsDeriv);
        eta0HzFourierP = eta0HzPrimaryFourier(n, K, W, x0, y0, z0, beta, bs);
        eta0HzFourierDerivP = eta0HzPrimaryFourierDeriv(n, K, W, x0, y0, z0, beta, bsDeriv);
        EphiFourierP = EphiPrimaryFourier(rho, n, K, W, EzFourierP, eta0HzFourierDerivP);
        ErhoFourierP = ErhoPrimaryFourier(rho, n, K, W, EzFourierDerivP, eta0HzFourierP);
        eta0HphiFourierP = eta0HphiPrimaryFourier(rho, n, K, W, EzFourierDerivP, eta0HzFourierP);
        eta0HrhoFourierP = eta0HrhoPrimaryFourier(rho, n, K, W, EzFourierP, eta0HzFourierDerivP);
    
        EzFourierS = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
        EphiFourierS = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        ErhoFourierS = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0HzFourierS = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0HphiFourierS = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0HrhoFourierS = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        
        EzFourier = EzFourierS - EzFourierP * heaviside(1 - rho);
        EphiFourier = EphiFourierS - EphiFourierP * heaviside(1 - rho);
        ErhoFourier = ErhoFourierS - ErhoFourierP * heaviside(1 - rho);
        eta0HzFourier = eta0HzFourierS - eta0HzFourierP * heaviside(1 - rho);
        eta0HphiFourier = eta0HphiFourierS - eta0HphiFourierP * heaviside(1 - rho);
        eta0HrhoFourier = eta0HrhoFourierS - eta0HrhoFourierP * heaviside(1 - rho);
           
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
