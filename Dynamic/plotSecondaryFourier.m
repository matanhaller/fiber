% Plotting Fourier transform of the secondary electric and magnetic fields
% and comparing to the static expression.

% Clearing workspace
close all;
clear;
clc;

epsilon = 12;
rho = 1.5;
x0 = 0; y0 = 2; z0 = 0;
omega = linspace(-6, 6, 3.2e4);
kz = 1;
beta = logspace(-2, log10(0.95), 10);
n = 1;

%% Intensity
figure; hold on;

for b=0.25
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, b, epsilon);
    EzFourier = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
    EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
    ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
    eta0HzFourier = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
    eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
    eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
    IE = abs(EzFourier).^2 + abs(EphiFourier).^2 + abs(ErhoFourier).^2;
    IH = abs(eta0HzFourier).^2 + abs(eta0HphiFourier).^2 + abs(eta0HrhoFourier).^2;
    
    plot(omega, abs(IH), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta = %.2f$', b));
end

xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Amplitude', 'FontSize', 14, 'Interpreter', 'latex');
title('Secondary Longitudinal Electric Field Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Integrating to check relation to delta function
figure; hold on;

for n=0:8
    disp(['n = ', int2str(n)]);
    intVec = zeros(1, numel(beta));
    intVecE = zeros(1, numel(beta));
    intVecH = zeros(1, numel(beta));
    for i=1:numel(beta)
        b = beta(i);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, b, epsilon);
        EzFourier = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
        EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
        ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
        eta0HzFourier = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
        eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
        eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
%         intVec(i) = trapz(omega, abs(EzFourier)).^2 + trapz(omega, abs(EphiFourier)).^2 + trapz(omega, abs(ErhoFourier)).^2 ...
%             + trapz(omega, abs(eta0HzFourier)).^2 + trapz(omega, abs(eta0HphiFourier)).^2 + trapz(omega, abs(eta0HrhoFourier)).^2;
%         intVecE(i) = trapz(omega, abs(EzFourier)).^2 + trapz(omega, abs(EphiFourier)).^2 + trapz(omega, abs(ErhoFourier)).^2;
%         intVecH(i) = trapz(omega, abs(eta0HzFourier)).^2 + trapz(omega, abs(eta0HphiFourier)).^2 + trapz(omega, abs(eta0HrhoFourier)).^2;
        intVecE(i) = abs(trapz(omega, EzFourier)).^2 + abs(trapz(omega, EphiFourier)).^2 + abs(trapz(omega, ErhoFourier)).^2;
        intVecH(i) = abs(trapz(omega, eta0HzFourier)).^2 + abs(trapz(omega, eta0HphiFourier)).^2 + abs(trapz(omega, eta0HrhoFourier)).^2;
    end
    
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
    
    Istatic = abs(EzFourierStatic).^2 + abs(EphiFourierStatic).^2 + abs(ErhoFourierStatic).^2;
    
    plot(beta, intVecH./intVecE, '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
       
end

xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$I_H/I_E$', 'FontSize', 14, 'Interpreter', 'latex');
title('Secondary Intensity, $\rho<R,k_z<\varepsilon\frac{\omega_\mathrm{max}^2}{c^2}$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Plotting field components

close all;
epsilon = 1.001;
omega = linspace(-10, 10, 100);
kz = 1;
beta = 0.8;
x0 = 0; y0 = 2; z0 = 1;
rho = 0.9;
n = 2;

[Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon);
EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, beta, epsilon);
EphiFourierP = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);

figure; hold on;
plot(omega, abs(EphiFourier), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
plot(omega, abs(EphiFourierP), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));

epsilon = 12;
n = 0;
kz = 0.001;
omega = 0.1;
rho = linspace(0.9, 1.1, 1000);
EzFourier = zeros(1, numel(rho));
EphiFourier = zeros(1, numel(rho));
ErhoFourier = zeros(1, numel(rho));
eta0HzFourier = zeros(1, numel(rho));
eta0HphiFourier = zeros(1, numel(rho));
eta0HrhoFourier = zeros(1, numel(rho));

[Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon);
for i=1:numel(rho)
    r = rho(i);
    EzFourier(i) = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, n, kz, omega, x0, y0, z0, beta, epsilon);
    EphiFourier(i) = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, n, kz, omega, x0, y0, z0, beta, epsilon);
    ErhoFourier(i) = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, n, kz, omega, x0, y0, z0, beta, epsilon);
    eta0HzFourier(i) = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, n, kz, omega, x0, y0, z0, beta, epsilon);
    eta0HphiFourier(i) = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, n, kz, omega, x0, y0, z0, beta, epsilon);
    eta0HrhoFourier(i) = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, n, kz, omega, x0, y0, z0, beta, epsilon);
end

EzFourierP = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
EphiFourierP = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
ErhoFourierP = ErhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
eta0HzFourierP = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
eta0HphiFourierP = eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
eta0HrhoFourierP = eta0HrhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);

figure; hold on;
plot(rho, abs(EzFourier + EzFourierP .* heaviside(rho - 1)), 'LineWidth', 1);

figure; hold on;
plot(rho, abs(EphiFourier + EphiFourierP .* heaviside(rho - 1)), 'LineWidth', 1);

figure; hold on;
plot(rho, abs(ErhoFourier + ErhoFourierP .* heaviside(rho - 1)), 'LineWidth', 1);

figure; hold on;
plot(rho, abs(eta0HzFourier + eta0HzFourierP .* heaviside(rho - 1)), 'LineWidth', 1);

figure; hold on;
plot(rho, abs(eta0HphiFourier + eta0HphiFourierP .* heaviside(rho - 1)), 'LineWidth', 1);

figure; hold on;
plot(rho, abs(eta0HrhoFourier + eta0HrhoFourierP .* heaviside(rho - 1)), 'LineWidth', 1);

% n = 1;
% kz = 0.01;
% omega = 0.00999;
% 
% drho = rho(2) - rho(1);
% Ez = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
% dEz = diff(Ez) / drho;
% dEz2 = EzPrimaryFourierDeriv(rho, n, kz, omega, x0, y0, z0, beta);
% 
% figure; hold on;
% plot(rho(1:end-1), abs(dEz), 'LineWidth', 1);
% plot(rho, abs(dEz2), 'LineWidth', 1);


% kz = 0.5;
% omega = 0.2;
% rho = 0.5;
% beta = logspace(-2, log10(0.95), 10);
% 
% figure; hold on;
% 
% for n=0:8
%     disp(['n = ', int2str(n)]);
%     EzFourier = zeros(1, numel(beta));
%     EphiFourier = zeros(1, numel(beta));
%     ErhoFourier = zeros(1, numel(beta));
%     eta0HzFourier = zeros(1, numel(beta));
%     eta0HphiFourier = zeros(1, numel(beta));
%     eta0HrhoFourier = zeros(1, numel(beta));
%     
%     for i=1:numel(beta)
%         b = beta(i);
%         [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, b, epsilon);
%         EzFourier(i) = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
%         EphiFourier(i) = 0*EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
%         ErhoFourier(i) = 0*ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
%         eta0HzFourier(i) = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
%         eta0HphiFourier(i) = 0*eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
%         eta0HrhoFourier(i) = 0*eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, b, epsilon);
%     end
%     
%     IE = abs(EzFourier).^2 + abs(EphiFourier).^2 + abs(ErhoFourier).^2;
%     IH = abs(eta0HzFourier).^2 + abs(eta0HphiFourier).^2 + abs(eta0HrhoFourier).^2;
%     
%     plot(beta, IH./IE, '--o', 'LineWidth', 1);
% end