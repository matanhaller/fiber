% Plotting Fourier transform of the primary electric and magnetic fields
% and comparing to the static expression.

% Clearing workspace
close all;
clear;
clc;

rho = 1;
x0 = 0; y0 = 1.05; z0 = 0;
kz = 0.5;
omega = logspace(log10(1.0001e-3), log10(20.0001), 8e3); omega = [-flip(omega), omega];
% omega = linspace(-10, 10, 200);
beta = linspace(1e-2, 0.2, 40);
n = 10;
nuMax = 100;
sigma = 0;

%% Plotting Fourier components
close all;

nuMax = 100;
sigma = 0;

M = 100;
kz = logspace(-3, log10(20), M); kz = [-flip(kz), kz];
omega = logspace(log10(1.0001e-3), log10(20.0001), M); omega = [-flip(omega), omega];

[K, W] = meshgrid(kz, omega);

tic;

beta = 0.99;

EzFourierSum = zeros(1, numel(omega));

for n=-20:20
    disp(n);
    bs = besselSum(rho, n, K, W, beta, nuMax, 0);
    bsDeriv = besselSumDeriv(rho, n, K, W, beta, nuMax, 0);    
    EzFourier = EzPrimaryFourier(n, K, W, x0, y0, z0, beta, bs, sigma);
    EzFourierSum = EzFourierSum + EzFourier;
end

figure;
surf(kz, omega, abs(EzFourierSum), 'EdgeColor', 'none');
view(2);
xlim([-20, 20]);
ylim([-20, 20]);
xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
colorbar;

toc;

%% Integrating to check relation to delta function of all fields
close all;

rho = 1;
x0 = 0; y0 = 1.05; z0 = 0;
kz = 10;
omega = logspace(log10(1.0001e-3), log10(20.0001), 8e3); omega = [-flip(omega), omega];
% omega = linspace(-10, 10, 200);
beta = linspace(1e-2, 0.2, 40);
nuMax = 100;
sigma = 0;

nuMax = 100;

tic;

sigma = 0;

for i=1:3
    figure(i); hold on;
end

for n=0:8
    fprintf('n = %d\n', n);
    
    EzFourierStatic = -(1/pi) * 1j*kz .* exp(1j*(pi/2)*n) .* besselk(n, abs(kz)*y0) .* besseli(n, abs(kz)*rho);
    EphiFourierStatic = -(1/pi) * (1j*n / rho) * exp(1j*(pi/2)*n) .* besselk(n, abs(kz)*y0) .* besseli(n, abs(kz)*rho);
    ErhoFourierStatic = (1/pi) * abs(kz) .* exp(1j*(pi/2)*n) .* besselk(n, abs(kz)*y0) .* besselip(n, abs(kz)*rho);

    intVecEz = zeros(1, numel(beta));
    intVecEphi = zeros(1, numel(beta));
    intVecErho = zeros(1, numel(beta));
    intVeceta0Hz = zeros(1, numel(beta));
    intVeceta0Hphi = zeros(1, numel(beta));
    intVeceta0Hrho = zeros(1, numel(beta));
    
    for i=1:numel(beta)
        b = beta(i);
        
        bs = besselSum(rho, n, kz, omega, b, nuMax, 0);
        bsDeriv = besselSumDeriv(rho, n, kz, omega, b, nuMax, 0);
        
        EzFourier = EzPrimaryFourier(n, kz, omega, x0, y0, z0, b, bs, sigma);
        EzFourierDeriv = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, b, bsDeriv, sigma);
        eta0HzFourier = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, b, bs, sigma);
        eta0HzFourierDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, b, bsDeriv, sigma);
        EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);
        ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
        eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
        eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);
        
        intVecEz(i) = trapz(omega, EzFourier);
        intVecEphi(i) = trapz(omega, EphiFourier);
        intVecErho(i) = trapz(omega, ErhoFourier);
        intVeceta0Hz(i) = trapz(omega, eta0HzFourier);
        intVeceta0Hphi(i) = trapz(omega, eta0HphiFourier);
        intVeceta0Hrho(i) = trapz(omega, eta0HrhoFourier);
    end
    
    figure(1);
    plot(beta, abs(intVecEz./EzFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(2);
    plot(beta, abs(intVecEphi./EphiFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(3);
    plot(beta, abs(intVecErho./ErhoFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
end

figure(1);
title('Longitudinal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(2);
title('Azimuthal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(3);
title('Radial Electric Field', 'FontSize', 14, 'Interpreter', 'latex');

for i=1:3
    figure(i);
    xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Ratio', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end

toc;

%% Plotting ratio between electric and magnetic field intensities
close all;

beta = 1 - 10.^(-linspace(1, 4, 5));
rho = 1;
x0 = 0; y0 = 1.05; z0 = 0;
kz = 1;
omega = logspace(log10(1.0001e-3), log10(20.0001), 8e3); omega = [-flip(omega), omega];
nuMax = 100;
n = 1;
sigma = 0;

tic;

for i=1:3
    figure(i); hold on;
end

for b=beta
    fprintf('beta = %d\n', b);
    
    bs = besselSum(rho, n, kz, omega, b, nuMax, 0);
    bsDeriv = besselSumDeriv(rho, n, kz, omega, b, nuMax, 0);
        
    EzFourier = EzPrimaryFourier(n, kz, omega, x0, y0, z0, b, bs, sigma);
    EzFourierDeriv = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, b, bsDeriv, sigma);
    eta0HzFourier = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, b, bs, sigma);
    eta0HzFourierDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, b, bsDeriv, sigma);
    EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);
    ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
    eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
    eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);
    
    figure(1);
    plot(omega, eta0HzFourier./EzFourier*(-1j)*sign(kz), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.5f$', b));
    figure(2);
    plot(omega, eta0HphiFourier./EphiFourier*(-1j)*sign(kz), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.5f$', b));
    figure(3);
    plot(omega, eta0HrhoFourier./ErhoFourier*(-1j)*sign(kz), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.5f$', b));
end
    
    

figure(1);
title('Longitudinal', 'FontSize', 14, 'Interpreter', 'latex');
figure(2);
title('Azimuthal', 'FontSize', 14, 'Interpreter', 'latex');
figure(3);
title('Radial', 'FontSize', 14, 'Interpreter', 'latex');

for i=1:3
    figure(i);
    xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Ratio', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end

toc;

%% Verifying Fourier transforms yield real functions
close all;

beta = 0.8;
n = -3;
rho = 1;
x0 = 0; y0 = 1.05; z0 = 0;
omega = logspace(log10(1.0001e-3), log10(20.0001), 8e3); omega = [-flip(omega), omega];
nuMax = 100;
sigma = 0;

kz = 10;

tic;

for i=1:6
    figure(i); hold on;
end

bs = besselSum(rho, n, kz, omega, beta, nuMax, 0);
bsDeriv = besselSumDeriv(rho, n, kz, omega, beta, nuMax, 0);

bs2 = besselSum(rho, -n, -kz, -omega, beta, nuMax, 0);
bsDeriv2 = besselSumDeriv(rho, -n, -kz, -omega, beta, nuMax, 0);
        
EzFourier = EzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
EzFourierDeriv = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
eta0HzFourier = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
eta0HzFourierDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);
ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);

EzFourier2 = EzPrimaryFourier(-n, -kz, -omega, x0, y0, z0, beta, bs2, sigma);
EzFourierDeriv2 = EzPrimaryFourierDeriv(-n, -kz, -omega, x0, y0, z0, beta, bsDeriv2, sigma);
eta0HzFourier2 = eta0HzPrimaryFourier(-n, -kz, -omega, x0, y0, z0, beta, bs2, sigma);
eta0HzFourierDeriv2 = eta0HzPrimaryFourierDeriv(-n, -kz, -omega, x0, y0, z0, beta, bsDeriv2, sigma);
EphiFourier2 = EphiPrimaryFourier(rho, -n, -kz, -omega, EzFourier2, eta0HzFourierDeriv2);
ErhoFourier2 = ErhoPrimaryFourier(rho, -n, -kz, -omega, EzFourierDeriv2, eta0HzFourier2);
eta0HphiFourier2 = eta0HphiPrimaryFourier(rho, -n, -kz, -omega, EzFourierDeriv2, eta0HzFourier2);
eta0HrhoFourier2 = eta0HrhoPrimaryFourier(rho, -n, -kz, -omega, EzFourier2, eta0HzFourierDeriv2);

EzSum = EzFourier + EzFourier2;
EphiSum = EphiFourier + EphiFourier2;
ErhoSum = ErhoFourier + ErhoFourier2;
eta0HzSum = eta0HzFourier + eta0HzFourier2;
eta0HphiSum = eta0HphiFourier + eta0HphiFourier2;
eta0HrhoSum = eta0HrhoFourier + eta0HrhoFourier2;

figure(1);
plot(omega, real(EzSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(EzSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(2);
plot(omega, real(EphiSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(EphiSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(3);
plot(omega, real(ErhoSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(ErhoSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(4);
plot(omega, real(eta0HzSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(eta0HzSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(5);
plot(omega, real(eta0HphiSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(eta0HphiSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(6);
plot(omega, real(eta0HrhoSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(eta0HrhoSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');

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
    xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Value', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end

toc;

%% Plotting inverse transform of the primary field
close all;

N = -10:10;
beta = 0.99;
x0 = 0; y0 = 1.4; z0 = 0;
nuMax = 100;
sigma = 0;

rhos = linspace(1e-3, 1.3, 40);
phi = linspace(-pi, pi, 241);

kz = linspace(-10, 10, 100);
omega = linspace(-10.0001, 10.0001, 100);
dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);
[K, W] = meshgrid(kz, omega);

[R, P] = meshgrid(rhos, phi);
X = R.*cos(P); Y = R.*sin(P);

Ez = zeros(size(X));
Ephi = zeros(size(X));
Erho = zeros(size(X));
eta0Hz = zeros(size(X));
eta0Hphi = zeros(size(X));
eta0Hrho = zeros(size(X));

tic;

for n=N
    disp(n);
    
    EzInv = zeros(1, numel(rhos));
    EphiInv = zeros(1, numel(rhos));
    ErhoInv = zeros(1, numel(rhos));
    eta0HzInv = zeros(1, numel(rhos));
    eta0HphiInv = zeros(1, numel(rhos));
    eta0HrhoInv = zeros(1, numel(rhos));
    
    for i=1:numel(rhos)
        rho = rhos(i);
        
        disp(rho);
        
        bs = besselSum(rho, n, K, W, beta, nuMax, 0);
        bsDeriv = besselSumDeriv(rho, n, K, W, beta, nuMax, 0);

        EzFourier = EzPrimaryFourier(n, K, W, x0, y0, z0, beta, bs, sigma);
        EzFourierDeriv = EzPrimaryFourierDeriv(n, K, W, x0, y0, z0, beta, bsDeriv, sigma);
        eta0HzFourier = eta0HzPrimaryFourier(n, K, W, x0, y0, z0, beta, bs, sigma);
        eta0HzFourierDeriv = eta0HzPrimaryFourierDeriv(n, K, W, x0, y0, z0, beta, bsDeriv, sigma);
        EphiFourier = EphiPrimaryFourier(rho, n, K, W, EzFourier, eta0HzFourierDeriv);
        ErhoFourier = ErhoPrimaryFourier(rho, n, K, W, EzFourierDeriv, eta0HzFourier);
        eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, K, W, EzFourierDeriv, eta0HzFourier);
        eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, K, W, EzFourier, eta0HzFourierDeriv);
        
        EzInv(i) = trapz(omega, trapz(kz, EzFourier, 2));
        EphiInv(i) = trapz(omega, trapz(kz, EphiFourier, 2));
        ErhoInv(i) = trapz(omega, trapz(kz, ErhoFourier, 2));
        eta0HzInv(i) = trapz(omega, trapz(kz, eta0HzFourier, 2));
        eta0HphiInv(i) = trapz(omega, trapz(kz, eta0HphiFourier, 2));
        eta0HrhoInv(i) = trapz(omega, trapz(kz, eta0HrhoFourier, 2));
    end
    
    Ez = Ez + exp(-1j*n*phi).' * EzInv;
    Ephi = Ephi + exp(-1j*n*phi).' * EphiInv;
    Erho = Erho + exp(-1j*n*phi).' * ErhoInv;
    eta0Hz = eta0Hz + exp(-1j*n*phi).' * eta0HzInv;
    eta0Hphi = eta0Hphi + exp(-1j*n*phi).' * eta0HphiInv;
    eta0Hrho = eta0Hrho + exp(-1j*n*phi).' * eta0HrhoInv;
end

for i=1:6
    figure(i); hold on;
end

figure(1);
surf(X, Y, real(Ez), 'EdgeColor', 'None');
view(2);
colorbar;

figure(2);
surf(X, Y, real(Ephi), 'EdgeColor', 'None');
view(2);
colorbar;

figure(3);
surf(X, Y, real(Erho), 'EdgeColor', 'None');
view(2);
colorbar;

figure(4);
surf(X, Y, real(eta0Hz), 'EdgeColor', 'None');
view(2);
colorbar;

figure(5);
surf(X, Y, real(eta0Hphi), 'EdgeColor', 'None');
view(2);
colorbar;

figure(6);
surf(X, Y, real(eta0Hrho), 'EdgeColor', 'None');
view(2);
colorbar;

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
    xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Value', 'FontSize', 14, 'Interpreter', 'latex');
end

toc;