% Plotting Fourier transform of the primary electric and magnetic fields
% and comparing to the static expression.

% Clearing workspace
close all;
clear;
clc;

rho = 0.5;
x0 = 0; y0 = 2; z0 = 0;
kz = 0.5;
omega = linspace(-4, 4, 1000);
beta = logspace(-3, -1, 10);
n = 0;

%% Longitudinal electric field
figure; hold on;

for b=beta
    primaryFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    plot(omega, abs(primaryFourier), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta = %.2f$', b));
end

xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Amplitude', 'FontSize', 14, 'Interpreter', 'latex');
title('Primary Longitudinal Electric Field Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

% Integrating to check relation to delta function
figure; hold on;

for n=0:5
    intVec = zeros(1, numel(beta));
    for i=1:numel(beta)
        b = beta(i);
        primaryFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        intVec(i) = trapz(omega, primaryFourier);
    end
    
    primaryFourierStatic = -(1/pi) * 1j*kz * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    
    plot(beta, abs(intVec)./abs(primaryFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
       
end

xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Relative Error', 'FontSize', 14, 'Interpreter', 'latex');
title('Comparison to Static Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Azimuthal electric field
figure; hold on;

for b=beta
    primaryFourier = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    plot(omega, abs(primaryFourier), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta = %.2f$', b));
end

xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Amplitude', 'FontSize', 14, 'Interpreter', 'latex');
title('Primary Azimuthal Electric Field Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

% Integrating to check relation to delta function
figure; hold on;

for n=0:5
    intVec = zeros(1, numel(beta));
    for i=1:numel(beta)
        b = beta(i);
        primaryFourier = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        intVec(i) = trapz(omega, primaryFourier);
    end
    
    primaryFourierStatic = -(1/pi) * (1j*n/rho) * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    
    plot(beta, abs(intVec)./abs(primaryFourierStatic), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
       
end

xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Relative Error', 'FontSize', 14, 'Interpreter', 'latex');
title('Comparison to Static Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Longitudinal magnetic field (times free-space wave impedance)
figure; hold on;

for b=beta
    primaryFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    plot(omega, abs(primaryFourier), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta = %.2f$', b));
end

xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Amplitude', 'FontSize', 14, 'Interpreter', 'latex');
title('Primary Longitudinal Magnetic Field Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

% Integrating to check relation to delta function
figure; hold on;

for n=0:5
    intVec = zeros(1, numel(beta));
    for i=1:numel(beta)
        b = beta(i);
        primaryFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        intVec(i) = trapz(omega, primaryFourier);
    end
    
    primaryFourierStatic = 0;
    
    plot(beta, abs(intVec), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
       
end

xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Relative Error', 'FontSize', 14, 'Interpreter', 'latex');
title('Comparison to Static Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Azimuthal magnetic field (times free-space wave impedance)
figure; hold on;

for b=beta
    primaryFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    plot(omega, abs(primaryFourier), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta = %.2f$', b));
end

xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Amplitude', 'FontSize', 14, 'Interpreter', 'latex');
title('Primary Azimuthal Magnetic Field Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

% Integrating to check relation to delta function
figure; hold on;

for n=0:5
    intVec = zeros(1, numel(beta));
    for i=1:numel(beta)
        b = beta(i);
        primaryFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        intVec(i) = trapz(omega, primaryFourier);
    end
    
    primaryFourierStatic = 0;
    
    plot(beta, abs(intVec), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
       
end

xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Relative Error', 'FontSize', 14, 'Interpreter', 'latex');
title('Comparison to Static Fourier Transform', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Ratio between longitudinal electric and magnetic field
n = 7;

figure; hold on;
for b=beta
    Ez = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    eta0Hz = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    plot(omega, abs(Ez./eta0Hz), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta = %.2f$', b));
end

xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Ratio', 'FontSize', 14, 'Interpreter', 'latex');
title('Ratio between Longitudinal Electric and Magnetic Field', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Ratio between azimuthal electric and magnetic field
n = 1;

figure; hold on;
for b=beta
    Ephi = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    eta0Hphi = eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
    plot(omega, abs(Ephi./eta0Hphi), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta = %.2f$', b));
end

xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Ratio', 'FontSize', 14, 'Interpreter', 'latex');
title('Ratio between Azimuthal Electric and Magnetic Field', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Intensity ratios
figure; hold on;

for n=0:8
    disp(['n = ', int2str(n)]);
    intVec = zeros(1, numel(beta));
    for i=1:numel(beta)
        b = beta(i);
        EzFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        eta0HzFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        intVec(i) = trapz(omega, abs(EzFourier)).^2 + trapz(omega, abs(EphiFourier)).^2 + trapz(omega, abs(ErhoFourier)).^2 ...
            + trapz(omega, abs(eta0HzFourier)).^2 + trapz(omega, abs(eta0HphiFourier)).^2 + trapz(omega, abs(eta0HrhoFourier)).^2;
    end
    
    EzFourierStatic = -(1/pi) * 1j*kz * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    EphiFourierStatic = -(1/pi) * (1j*n / rho) * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    ErhoFourierStatic = (1/pi) * 1j*abs(kz) * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besselip(n, abs(kz)*rho);
    Istatic = abs(EzFourierStatic).^2 + abs(EphiFourierStatic).^2 + abs(ErhoFourierStatic).^2;

    plot(beta, intVec ./ Istatic, '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
       
end

xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$I_\mathrm{dynamic}/I_\mathrm{static}$', 'FontSize', 14, 'Interpreter', 'latex');
title('Comparison to Static Fourier Transform - Primary Intensity', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Ratio beween electric and magnetic field intensities
beta = logspace(-2, log10(0.999), 20);

figure; hold on;

for n=0:8
    disp(['n = ', int2str(n)]);
    intVecE = zeros(1, numel(beta));
    intVecH = zeros(1, numel(beta));
    for i=1:numel(beta)
        b = beta(i);
        EzFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        eta0HzFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, b);
        intVecE(i) = trapz(omega, abs(EzFourier)).^2 + trapz(omega, abs(EphiFourier)).^2 + trapz(omega, abs(ErhoFourier)).^2;
        intVecH(i) = trapz(omega, abs(eta0HzFourier)).^2 + trapz(omega, abs(eta0HphiFourier)).^2 + trapz(omega, abs(eta0HrhoFourier)).^2;
    end
    
    EzFourierStatic = -(1/pi) * 1j*kz * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    EphiFourierStatic = -(1/pi) * (1j*n / rho) * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besseli(n, abs(kz)*rho);
    ErhoFourierStatic = (1/pi) * 1j*abs(kz) * exp(1j*(pi/2)*n) * besselk(n, abs(kz)*y0) * besselip(n, abs(kz)*rho);
    Istatic = abs(EzFourierStatic).^2 + abs(EphiFourierStatic).^2 + abs(ErhoFourierStatic).^2;

    plot(beta, intVecH ./ intVecE, '--o', 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
       
end

xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$I_H/I_E$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');