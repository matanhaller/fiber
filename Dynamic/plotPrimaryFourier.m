% Plotting Fourier transform of the primary electric and magnetic fields
% and comparing to the static expression.

% Clearing workspace
close all;
clear;
clc;

rho = 1.25;
x0 = 0; y0 = 2; z0 = 0;
kz = 2;
omega = linspace(-10, 10, 8e3);
beta = linspace(1e-2, 0.95, 10);
n = -5;
nuMax = 40;

%% Integrating to check relation to delta function of all fields
tic;

for i=1:4
    figure(i); hold on;
end

for n=0:5
    fprintf('n = %d\n', n);
    
    EzFourierStatic = -(1/pi) * 1j*kz * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    EphiFourierStatic = -(1/pi) * (1j*n / rho) * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    ErhoFourierStatic = (1/pi) * 1j*abs(kz) * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besselip(n, abs(kz)*rho);
    
    intVecEz = zeros(1, numel(beta));
    intVecEphi = zeros(1, numel(beta));
    intVecErho = zeros(1, numel(beta));
    intVeceta0Hz = zeros(1, numel(beta));
    intVeceta0Hphi = zeros(1, numel(beta));
    intVeceta0Hrho = zeros(1, numel(beta));
    
    for i=1:numel(beta)
        b = beta(i);
        
        bs = besselSum(rho, n, kz, omega, b, nuMax);
        bsDeriv = besselSumDeriv(rho, n, kz, omega, b, nuMax);
        
        EzFourier = EzPrimaryFourier(n, kz, omega, x0, y0, z0, b, bs);
        EzFourierDeriv = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, b, bsDeriv);
        eta0HzFourier = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, b, bs);
        eta0HzFourierDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, b, bsDeriv);
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
    plot(beta, abs(intVecEz)./abs(EzFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(2);
    plot(beta, abs(intVecEphi)./abs(EphiFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(3);
    plot(beta, abs(intVecErho)./abs(ErhoFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(4);
    plot(beta, (abs(intVeceta0Hz).^2 + abs(intVeceta0Hphi).^2 + abs(intVeceta0Hrho).^2) ...
        ./ (abs(intVecEz).^2 + abs(intVecEphi).^2 + abs(intVecErho).^2), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
end

figure(1);
title('Longitudinal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(2);
title('Azimuthal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(3);
title('Radial Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(4);
title('Ratio between $E$ and $H$ Intensities', 'FontSize', 14, 'Interpreter', 'latex');

for i=1:4
    figure(i);
    xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Ratio', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end

toc;

%% Verifying Fourier transforms yield real functions
tic;

for i=1:6
    figure(i); hold on;
end

beta = 0.8;
n = 8;
kz = 0.01;

bs = besselSum(rho, n, kz, omega, beta, nuMax);
bsDeriv = besselSumDeriv(rho, n, kz, omega, beta, nuMax);

bs2 = besselSum(rho, -n, -kz, -omega, beta, nuMax);
bsDeriv2 = besselSumDeriv(rho, -n, -kz, -omega, beta, nuMax);
        
EzFourier = EzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs);
EzFourierDeriv = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv);
eta0HzFourier = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs);
eta0HzFourierDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv);
EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);
ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, EzFourierDeriv, eta0HzFourier);
eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, EzFourier, eta0HzFourierDeriv);

EzFourier2 = EzPrimaryFourier(-n, -kz, -omega, x0, y0, z0, beta, bs2);
EzFourierDeriv2 = EzPrimaryFourierDeriv(-n, -kz, -omega, x0, y0, z0, beta, bsDeriv2);
eta0HzFourier2 = eta0HzPrimaryFourier(-n, -kz, -omega, x0, y0, z0, beta, bs2);
eta0HzFourierDeriv2 = eta0HzPrimaryFourierDeriv(-n, -kz, -omega, x0, y0, z0, beta, bsDeriv2);
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
