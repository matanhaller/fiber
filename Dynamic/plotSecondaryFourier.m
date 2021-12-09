% Plotting Fourier transform of the secondary electric and magnetic fields
% and comparing to the static expression.

% Clearing workspace
close all;
clear;
clc;

epsilon = 1.2;
rho = 1.2;
x0 = 0; y0 = 1.4; z0 = 0;
omega = linspace(-10, 10, 8e3);
kz = 0.1;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = (0.1:0.1:1) * gbC; 
% beta = sqrt(gb.^2 ./ (gb.^2 + 1));
beta = linspace(1e-2, 0.2, 40);
% beta = logspace(-2, log10(0.95), 20);
n = -7;
nuMax = 100;

%% Plotting Fourier components
% close all;

sigma = 0.0;
nuMax = 100;
epsilon = 12;
eps = 1e-1;

y0 = 1.05;
rho = 0.1;

Mw = 4e2; Mk = 4e2;

omegaMax = 10;
kzMax = 20;
omega = linspace(-omegaMax, omegaMax, Mw);
kz = linspace(-kzMax, kzMax, Mk);
[W, K] = meshgrid(omega, kz);

tic;

beta = 0.2;

AnkSum = zeros((size(W)));
BnkSum = zeros((size(W)));
eta0CnkSum = zeros((size(W)));
eta0DnkSum = zeros((size(W)));
EzFourierSum = zeros(size(W));
EphiFourierSum = zeros(size(W));
ErhoFourierSum = zeros(size(W));
eta0HzFourierSum = zeros(size(W));
eta0HphiFourierSum = zeros(size(W));
eta0HrhoFourierSum = zeros(size(W));

for n=0
    disp(n);
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0, eps);
   
    EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon, eps);
    EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);
    ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);
    eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);
    eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);
    eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);

%     eta0HzFourier(filterIdx) = 0;
    
    AnkSum = AnkSum + Ank;
    BnkSum = BnkSum + Bnk;
    eta0CnkSum = eta0CnkSum + eta0Cnk;
    eta0DnkSum = eta0DnkSum + eta0Dnk;

    EzFourierSum = EzFourierSum + EzFourier;
    EphiFourierSum = EphiFourierSum + EphiFourier;
    ErhoFourierSum = ErhoFourierSum + ErhoFourier;
    eta0HzFourierSum = eta0HzFourierSum + eta0HzFourier;
    eta0HphiFourierSum = eta0HphiFourierSum + eta0HphiFourier;
    eta0HrhoFourierSum = eta0HrhoFourierSum + eta0HrhoFourier;
end

figure; hold on;
surf(W, K, abs(eta0HzFourierSum), 'EdgeColor', 'none');
view(2);
xlim([omega(1), omega(end)]);
ylim([kz(1), kz(end)]);
xlabel('$\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{k}_z$', 'FontSize', 14, 'Interpreter', 'latex');
colorbar;
% caxis([0, 0.15]);

x = linspace(-20, 20, 1e3);
plot3(x, x, 0.1 * ones(1, numel(x)), '--w', 'LineWidth', 1);
plot3(x, -x, 0.1 * ones(1, numel(x)), '--w', 'LineWidth', 1);
plot3(x, x * sqrt(epsilon), 0.1 * ones(1, numel(x)), '--w', 'LineWidth', 1);
plot3(x, -x * sqrt(epsilon), 0.1 * ones(1, numel(x)), '--w', 'LineWidth', 1);

toc;

%% Integrating to check relation to delta function of all fields
tic;

esilon = 12;
sigma = 0;
eps = 0;

for i=1:3
    figure(i); hold on;
end

for n=0:5
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
        
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, b, epsilon, nuMax, sigma, 0, eps);
        
        EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, kz, omega, epsilon, eps);
        EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps);
        ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps);
        eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps);
        eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps);
        eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps);

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

%% Plotting ratio between electric and magnetic fields in the far field
close all;


sigma = 0;
nuMax = 100;
epsilon = 4;

rhos = logspace(log10(1.02), log10(50), 1e3);
beta = 0.1;

omega = 1;
kz = 0.95 * omega;

EzFourier = zeros(1, numel(rhos));
EphiFourier = zeros(1, numel(rhos));
ErhoFourier = zeros(1, numel(rhos));
eta0HzFourier = zeros(1, numel(rhos));
eta0HphiFourier = zeros(1, numel(rhos));
eta0HrhoFourier = zeros(1, numel(rhos));

figure; hold on;
for n=0:8
    for i=1:numel(rhos)
        rho = rhos(i);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);
                
        EzFourier(i) = EzSecondaryFourier(Ank, Bnk, rho, n, kz, omega, epsilon);
        EphiFourier(i) = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        ErhoFourier(i) = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        eta0HzFourier(i) = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        eta0HphiFourier(i) = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
        eta0HrhoFourier(i) = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    end
    
    EFourier = abs(EzFourier).^2 + abs(EphiFourier).^2 + abs(ErhoFourier).^2;
    eta0HFourier = abs(eta0HzFourier).^2 + abs(eta0HphiFourier).^2 + abs(eta0HrhoFourier).^2;

    plot(rhos, eta0HFourier./EFourier, 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
end

xlabel('$\rho$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Ratio', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

% set(gca, 'XScale', 'log');

%% Verifying Fourier transforms yield real functions
close all;

epsilon = 12;
omega = linspace(2e-2, 12, 100);
beta = 0.99;
x0 = 0; y0 = 1.4; z0 = 0;
rho = 1.05;
n = 10;
kz = 0.1;
nuMax = 100;
sigma = 0;

tic;

for i=1:6
    figure(i); hold on;
end

[Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);
[Ank2, Bnk2, eta0Cnk2, eta0Dnk2] = secondaryFieldCoeffs(-n, -kz, -omega, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);

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
omega = 2;
beta = 0.8;
x0 = 0; y0 = 1.4; z0 = 0;
n = -1;
kz = 0.2;
nuMax = 100;
sigma = 0.0;

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
    
    bs = besselSum(rho, n, kz, omega, beta, nuMax, 0);
    bsDeriv = besselSumDeriv(rho, n, kz, omega, beta, nuMax, 0);

    EzFourierP = EzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
    EzFourierDerivP = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
    eta0HzFourierP = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
    eta0HzFourierDerivP = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
    EphiFourierP = EphiPrimaryFourier(rho, n, kz, omega, EzFourierP, eta0HzFourierDerivP);
    ErhoFourierP = ErhoPrimaryFourier(rho, n, kz, omega, EzFourierDerivP, eta0HzFourierP);
    eta0HphiFourierP = eta0HphiPrimaryFourier(rho, n, kz, omega, EzFourierDerivP, eta0HzFourierP);
    eta0HrhoFourierP = eta0HrhoPrimaryFourier(rho, n, kz, omega, EzFourierP, eta0HzFourierDerivP);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax, sigma, 0, Inf);

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

%% Verifying Fourier transforms are symmetric w.r.t. kz
close all;

epsilon = 12;
omega = linspace(-10, 10, 200);
beta = 0.99;
x0 = 0; y0 = 1.4; z0 = 0;
rho = 0.25;
n = -1;
kz = 1e-3;
nuMax = 100;
sigma = 0;

tic;

for i=1:6
    figure(i); hold on;
end

[Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);
[Ank2, Bnk2, eta0Cnk2, eta0Dnk2] = secondaryFieldCoeffs(n, -kz, omega, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);

EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, kz, omega, epsilon);
EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);

EzFourier2 = EzSecondaryFourier(Ank2, Bnk2, rho, n, -kz, omega, epsilon);
EphiFourier2 = EphiSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, n, -kz, omega, epsilon);
ErhoFourier2 = ErhoSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, n, -kz, omega, epsilon);
eta0HzFourier2 = eta0HzSecondaryFourier(eta0Cnk2, eta0Dnk2, rho, n, -kz, omega, epsilon);
eta0HphiFourier2 = eta0HphiSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, n, -kz, omega, epsilon);
eta0HrhoFourier2 = eta0HrhoSecondaryFourier(Ank2, Bnk2, eta0Cnk2, eta0Dnk2, rho, n, -kz, omega, epsilon);

EzFourierSum = 1 + EzFourier2./EzFourier; % Anti-symmetric
EphiFourierSum = 1 - EphiFourier2./EphiFourier; % Symmetric
ErhoFourierSum = 1 - ErhoFourier2./ErhoFourier; % Symmetric
eta0HzFourierSum = 1 - eta0HzFourier2./eta0HzFourier; % Symmetric
eta0HphiFourierSum = 1 + eta0HphiFourier2./eta0HphiFourier; % Anti-symmetric
eta0HrhoFourierSum = 1 + eta0HrhoFourier2./eta0HrhoFourier; % Anti-symmetric

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

%% Plotting contribution of each harmonic
close all;

epsilon = 4;
beta = 0.8;
x0 = 0; y0 = 1.4; z0 = 0;
rho = 0.75;
N = 0:20;
M = 800;
kz = linspace(-10, 10, M);
dkz = kz(2) - kz(1);
omega = linspace(-10.0001, 10.0001, M);
domega = omega(2) - omega(1);
[K, W] = meshgrid(kz, omega);
nuMax = 100;
sigma = 0;

tic;

for i=1:6
    figure(i); hold on;
end

EzNorm = zeros(1, numel(N));
EphiNorm = zeros(1, numel(N));
ErhoNorm = zeros(1, numel(N));
eta0HzNorm = zeros(1, numel(N));
eta0HphiNorm = zeros(1, numel(N));
eta0HrhoNorm = zeros(1, numel(N));

parfor i=1:numel(N)
    n = N(i);
    disp(n);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);

    Ez = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
    Ez = sum(sum(Ez)) * dkz * domega;
    Ephi = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
    Ephi = sum(sum(Ephi)) * dkz * domega;
    Erho = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
    Erho = sum(sum(Erho)) * dkz * domega;
    eta0Hz = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
    eta0Hz = sum(sum(eta0Hz)) * dkz * domega;
    eta0Hphi = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
    eta0Hphi = sum(sum(eta0Hphi)) * dkz * domega;
    eta0Hrho = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
    eta0Hrho = sum(sum(eta0Hrho)) * dkz * domega;

    EzNorm(i) = norm(Ez, 'fro');
    EphiNorm(i) = norm(Ephi, 'fro');
    ErhoNorm(i) = norm(Erho, 'fro');
    eta0HzNorm(i) = norm(eta0Hz, 'fro');
    eta0HphiNorm(i) = norm(eta0Hphi, 'fro');
    eta0HrhoNorm(i) = norm(eta0Hrho, 'fro');
end

figure(1);
stem(N, EzNorm, 'LineWidth', 1);
figure(2);
stem(N, EphiNorm, 'LineWidth', 1);
figure(3);
stem(N, ErhoNorm, 'LineWidth', 1);
figure(4);
stem(N, eta0HzNorm, 'LineWidth', 1);
figure(5);
stem(N, eta0HphiNorm, 'LineWidth', 1);
figure(6);
stem(N, eta0HrhoNorm, 'LineWidth', 1);


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
    xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Value', 'FontSize', 14, 'Interpreter', 'latex');
end

toc;

%% Verifying kz and omega range and resolution are sufficient
close all;

M = 200;
epsilon = 12;
kz = linspace(-40, 40, M);
omega = linspace(-10.001, 10.001, M);
[K, W] = meshgrid(kz, omega);
beta = 0.99;
x0 = 0; y0 = 1.05; z0 = 0;
rho = 1.25;
n = -7;
nuMax = 40;
sigma = 0;

tic;

for i=1:6
    figure(i); hold on;
end

[Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);

EzFourier = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);

figure(1);
plot(kz, abs(EzFourier(0.5*numel(omega)+1,:)), 'LineWidth', 1);
figure(2);
plot(kz, abs(EphiFourier(0.5*numel(omega)+1,:)), 'LineWidth', 1);
figure(3);
plot(kz, abs(ErhoFourier(0.5*numel(omega)+1,:)), 'LineWidth', 1);
figure(4);
plot(kz, abs(eta0HzFourier(0.5*numel(omega)+1,:)), 'LineWidth', 1);
figure(5);
plot(kz, abs(eta0HphiFourier(0.5*numel(omega)+1,:)), 'LineWidth', 1);
figure(6);
plot(kz, abs(eta0HrhoFourier(0.5*numel(omega)+1,:)), 'LineWidth', 1);


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
    xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Value', 'FontSize', 14, 'Interpreter', 'latex');
end

toc;
