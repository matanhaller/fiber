% Plotting Fourier transform of the secondary electric and magnetic fields
% and comparing to the static expression.

% Clearing workspace
close all;
clear;
clc;

epsilon = 2;
rho = 0.5;
x0 = 0; y0 = 1.4; z0 = 0;
omega = linspace(-12, 12, 1e3);
kz = 1;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = (0.1:0.1:2) * gbC; 
beta = sqrt(gb.^2 ./ (gb.^2 + 1));
% beta = linspace(1e-2, 0.99, 40);
% beta = logspace(-2, log10(0.95), 20);
n = -3;
nuMax = 200;

%% Integrating to check relation to delta function of all fields
tic;

for i=1:4
    figure(i); hold on;
end

for n=1
    fprintf('n = %d\n', n);
    
    intVecEz = zeros(1, numel(beta));
    intVecEphi = zeros(1, numel(beta));
    intVecErho = zeros(1, numel(beta));
    intVeceta0Hz = zeros(1, numel(beta));
    intVeceta0Hphi = zeros(1, numel(beta));
    intVeceta0Hrho = zeros(1, numel(beta));
    
    BnkStatic = -(epsilon - 1) * (besselk(n,abs(kz)*y0) .* besselip(n,abs(kz)) .* besseli(n,abs(kz))) ...
        ./ (epsilon .* besselip(n,abs(kz)) .* besselk(n,abs(kz)) - besseli(n,abs(kz)) .* besselkp(n,abs(kz)));

    AnkStatic = (epsilon./besseli(n,abs(kz))) .* (besselk(n,abs(kz)*abs(y0)).*besseli(n,abs(kz)) + BnkStatic.*besselk(n,abs(kz)));

    PotentialSecondaryStatic = -(1/pi) * exp(1j*(pi/2)*n) ...
        * ((1/epsilon) .* AnkStatic .* besseli(n,abs(kz)*rho) * heaviside(1 - rho) ...
        + BnkStatic .* besselk(n,abs(kz)*rho) * heaviside(rho - 1));
    
    EzFourierStatic = - PotentialSecondaryStatic .* (-1j*kz);
    EphiFourierStatic = - (1./rho) .* PotentialSecondaryStatic .* (-1j*n);
    ErhoFourierStatic = (1/pi) * exp(1j*(pi/2)*n) ...
        * ((1/epsilon) .* AnkStatic .* abs(kz) .* besselip(n,abs(kz)*rho) * heaviside(1 - rho) ...
        + BnkStatic .* abs(kz) .* besselkp(n,abs(kz)*rho) * heaviside(rho - 1));
    
    for i=1:numel(beta)
        b = beta(i);
        
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, b, epsilon, nuMax);
        
        EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, kz, omega, epsilon);
        EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);

        intVecEz(i) = abs(trapz(omega, EzFourier));
        intVecEphi(i) = abs(trapz(omega, EphiFourier));
        intVecErho(i) = abs(trapz(omega, ErhoFourier));
        intVeceta0Hz(i) = abs(trapz(omega, eta0HzFourier));
        intVeceta0Hphi(i) = abs(trapz(omega, eta0HphiFourier));
        intVeceta0Hrho(i) = abs(trapz(omega, eta0HrhoFourier));
    end
    
    figure(1);
    plot(beta, abs(intVecEz)./abs(EzFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(2);
    plot(beta, abs(intVecEphi)./abs(EphiFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(3);
    plot(beta, abs(intVecErho)./abs(ErhoFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    figure(4);
    plot(beta, (abs(intVeceta0Hz).^2 + abs(intVeceta0Hphi).^2 + abs(intVeceta0Hrho).^2) ...
        ./(abs(intVecEz).^2 + abs(intVecEphi).^2 + abs(intVecErho).^2), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));   
end

figure(1);
title('Longitudinal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(2);
title('Azimuthal Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(3);
title('Radial Electric Field', 'FontSize', 14, 'Interpreter', 'latex');
figure(4);
title('Ratio between $E$ and $H$ intensities', 'FontSize', 14, 'Interpreter', 'latex');

for i=1:4
    figure(i);
    xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Ratio', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end

toc;

%% Verifying Fourier transforms yield real functions
close all;

epsilon = 1.2;
omega = linspace(-10, 10, 200);
beta = 0.8;
x0 = 0; y0 = 1.4; z0 = 0;
rho = 1.25;
n = -6;
kz = 3;
nuMax = 40;

tic;

for i=1:6
    figure(i); hold on;
end

[Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax);
[Ank2, Bnk2, eta0Cnk2, eta0Dnk2] = secondaryFieldCoeffs(-n, -kz, -omega, x0, y0, z0, beta, epsilon, nuMax);

EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, kz, omega, epsilon);
EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);

EzFourier2 = EzSecondaryFourier(Ank2, Bnk2, rho, -n, -kz, -omega, epsilon);
EphiFourier2 = EphiSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, -n, -kz, -omega, epsilon);
ErhoFourier2 = ErhoSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, -n, -kz, -omega, epsilon);
eta0HzFourier2 = eta0HzSecondaryFourier(eta0Cnk2, eta0Dnk2, rho, -n, -kz, -omega, epsilon);
eta0HphiFourier2 = eta0HphiSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, -n, -kz, -omega, epsilon);
eta0HrhoFourier2 = eta0HrhoSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, -n, -kz, -omega, epsilon);

EzFourierSum = EzFourier + EzFourier2;
EphiFourierSum = EphiFourier + EphiFourier2;
ErhoFourierSum = ErhoFourier + ErhoFourier2;
eta0HzFourierSum = eta0HzFourier + eta0HzFourier2;
eta0HphiFourierSum = eta0HphiFourier + eta0HphiFourier2;
eta0HrhoFourierSum = eta0HrhoFourier + eta0HrhoFourier2;

figure(1);
plot(omega, real(EzFourierSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(EzFourierSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(2);
plot(omega, real(EphiFourierSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(EphiFourierSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(3);
plot(omega, real(ErhoFourierSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(ErhoFourierSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(4);
plot(omega, real(eta0HzFourierSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(eta0HzFourierSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(5);
plot(omega, real(eta0HphiFourierSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(eta0HphiFourierSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');
figure(6);
plot(omega, real(eta0HrhoFourierSum), 'LineWidth', 1, 'DisplayName', 'Real');
plot(omega, imag(eta0HrhoFourierSum), 'LineWidth', 1, 'DisplayName', 'Imaginary');

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

%% Verifying boundary conditions

close all;
epsilon = 12;
omega = 9;
beta = 0.8;
x0 = 0; y0 = 1.4; z0 = 0;
n = -3;
kz = 0.1;
nuMax = 40;

rhos = linspace(0.5, 1.39, 1e3);

tic;

for i=1:6
    figure(i); hold on;
end

EzFourier = zeros(1, numel(rhos));
EphiFourier = zeros(1, numel(rhos));
ErhoFourier = zeros(1, numel(rhos));
eta0HzFourier = zeros(1, numel(rhos));
eta0HphiFourier = zeros(1, numel(rhos));
eta0HrhoFourier = zeros(1, numel(rhos));

for i=1:numel(rhos)
    rho = rhos(i);
    
    bs = besselSum(rho, n, kz, omega, beta, nuMax);
    bsDeriv = besselSumDeriv(rho, n, kz, omega, beta, nuMax);

    EzFourierP = EzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs);
    EzFourierDerivP = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv);
    eta0HzFourierP = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs);
    eta0HzFourierDerivP = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv);
    EphiFourierP = EphiPrimaryFourier(rho, n, kz, omega, EzFourierP, eta0HzFourierDerivP);
    ErhoFourierP = ErhoPrimaryFourier(rho, n, kz, omega, EzFourierDerivP, eta0HzFourierP);
    eta0HphiFourierP = eta0HphiPrimaryFourier(rho, n, kz, omega, EzFourierDerivP, eta0HzFourierP);
    eta0HrhoFourierP = eta0HrhoPrimaryFourier(rho, n, kz, omega, EzFourierP, eta0HzFourierDerivP);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax);

    EzFourierS = EzSecondaryFourier(Ank, Bnk, rho, n, kz, omega, epsilon);
    EphiFourierS = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    ErhoFourierS = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    eta0HzFourierS = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    eta0HphiFourierS = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    eta0HrhoFourierS = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    
    EzFourier(i) = EzFourierS + EzFourierP * heaviside(rho - 1);
    EphiFourier(i) = EphiFourierS + EphiFourierP * heaviside(rho - 1);
    ErhoFourier(i) = ErhoFourierS * epsilon * heaviside(1 - rho) + (ErhoFourierS + ErhoFourierP) * heaviside(rho - 1);
    eta0HzFourier(i) = eta0HzFourierS + eta0HzFourierP * heaviside(rho - 1);
    eta0HphiFourier(i) = eta0HphiFourierS + eta0HphiFourierP * heaviside(rho - 1);
    eta0HrhoFourier(i) = eta0HrhoFourierS + eta0HrhoFourierP * heaviside(rho - 1);
end

figure(1);
plot(rhos, abs(EzFourier), 'LineWidth', 1);
figure(2);
plot(rhos, abs(EphiFourier), 'LineWidth', 1);
figure(3);
plot(rhos, abs(ErhoFourier), 'LineWidth', 1);
figure(4);
plot(rhos, abs(eta0HzFourier), 'LineWidth', 1);
figure(5);
plot(rhos, abs(eta0HphiFourier), 'LineWidth', 1);
figure(6);
plot(rhos, abs(eta0HrhoFourier), 'LineWidth', 1);

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
    xlabel('$\rho$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Value', 'FontSize', 14, 'Interpreter', 'latex');
end

toc;
