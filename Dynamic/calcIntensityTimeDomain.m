% Calculating intensity in the time domain.

% Clearing workspace
close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 0.99;

sigma = 0;
nuMax = 100;

M = 120;
rhos = linspace(1e-3, 5, M);

Nk = 400;
Nw = 400;
kz = linspace(-10.0001, 10.0001, Nk);
omega = linspace(-20, 20, Nw);
[K, W] = meshgrid(kz, omega);

dkz = kz(2) - kz(1);
dw = omega(2) - omega(1);

z = linspace(-pi/dkz, pi/dkz, Nk+1); z(end) = [];
t = linspace(-pi/dw, pi/dw, Nw+1); t(1) = [];

tic
for n=-20:20
    disp(n);
    
    EzTotal = zeros(M, Nw);
    EphiTotal = zeros(M, Nw);
    ErhoTotal = zeros(M, Nw);
    eta0HzTotal = zeros(M, Nw);
    eta0HphiTotal = zeros(M, Nw);
    eta0HrhoTotal = zeros(M, Nw);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0, 1e2);
    
    for i=1:numel(rhos)
        rho = rhos(i);
        
        bs = besselSum(rho, n, K, W, beta, nuMax, 0);
        bsDeriv = besselSumDeriv(rho, n, K, W, beta, nuMax, 0);

        EzFourierP = EzPrimaryFourier(n, K, W, x0, y0, z0, beta, bs, sigma);
        EzFourierDerivP = EzPrimaryFourierDeriv(n, K, W, x0, y0, z0, beta, bsDeriv, sigma);
        eta0HzFourierP = eta0HzPrimaryFourier(n, K, W, x0, y0, z0, beta, bs, sigma);
        eta0HzFourierDerivP = eta0HzPrimaryFourierDeriv(n, K, W, x0, y0, z0, beta, bsDeriv, sigma);
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

        EzFourier = EzFourier .* (abs(W) - abs(K) > 0);
        EphiFourier = EphiFourier .* (abs(W) - abs(K) > 0);
        ErhoFourier = ErhoFourier .* (abs(W) - abs(K) > 0);
        eta0HzFourier = eta0HzFourier .* (abs(W) - abs(K) > 0);
        eta0HphiFourier = eta0HphiFourier .* (abs(W) - abs(K) > 0);
        eta0HrhoFourier = eta0HrhoFourier .* (abs(W) - abs(K) > 0);
        
        EzFFT = trapz(dkz, EzFourier, 2);
        EphiFFT = trapz(dkz, EphiFourier, 2);
        ErhoFFT = trapz(dkz, ErhoFourier, 2);
        eta0HzFFT = trapz(dkz, eta0HzFourier, 2);
        eta0HphiFFT = trapz(dkz, eta0HphiFourier, 2);
        eta0HrhoFFT = trapz(dkz, eta0HrhoFourier, 2);

        EzIFFT = fftshift(fft(fftshift(EzFFT))) * dw;
        EphiIFFT = fftshift(fft(fftshift(EphiFFT))) * dw;
        ErhoIFFT = fftshift(fft(fftshift(ErhoFFT))) * dw;
        eta0HzIFFT = fftshift(fft(fftshift(eta0HzFFT))) * dw;
        eta0HphiIFFT = fftshift(fft(fftshift(eta0HphiFFT))) * dw;
        eta0HrhoIFFT = fftshift(fft(fftshift(eta0HrhoFFT))) * dw;

        EzTotal(i,:) = single(flip(EzIFFT)).';
        EphiTotal(i,:) = single(flip(EphiIFFT)).';
        ErhoTotal(i,:) = single(flip(ErhoIFFT)).';
        eta0HzTotal(i,:) = single(flip(eta0HzIFFT)).';
        eta0HphiTotal(i,:) = single(flip(eta0HphiIFFT)).';
        eta0HrhoTotal(i,:) = single(flip(eta0HrhoIFFT)).';
        
        disp(rho);

    end
    save(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Ez_n=%d.mat', n), 'EzTotal');
    save(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Ephi_n=%d.mat', n), 'EphiTotal');
    save(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Erho_n=%d.mat', n), 'ErhoTotal');
    save(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hz_n=%d.mat', n), 'eta0HzTotal');
    save(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hphi_n=%d.mat', n), 'eta0HphiTotal');
    save(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hrho_n=%d.mat', n), 'eta0HrhoTotal');
end
toc
