% Plotting determinant of the secondary fields' coefficient matrix.

% Clearing workspace
close all;
clear;
clc;

epsilon = 12; % Relative perimittivity
M = 1e4;
n = 0;

%% Plotting for kz=0 (Cutoff frequency)
kz = 0;
omega = linspace(0, 4.5, M+1);

% Dispersion relation wavevectors
kVac = sqrt(kz.^2 - omega.^2);    
kCyl = sqrt(kz.^2 - epsilon * omega.^2);
       
% Determinant
I = besseli(n,kCyl);
Ip = besselip(n,kCyl);
K = besselk(n,kVac);
Kp = besselkp(n,kVac);
Delta = omega.^2 .* ((1./kVac).*(Kp./K) - epsilon.*(1./kCyl).*(Ip./I)) ...
    .* ((1./kVac).*(Kp./K) - (1./kCyl).*(Ip./I)) ...
    - n^2.*kz.^2.*(1./(kVac.^2) - 1./(kCyl.^2));

figure; hold on;
plot(omega, real(Delta), 'LineWidth', 1);
plot(omega, imag(Delta), 'LineWidth', 1);
ylim([-10, 10]);
xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\Delta/K_n$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting for specific omega
close all;

n = 2;
omega = 1.4;
kz = linspace(0, 10, M+1);

% Dispersion relation wavevectors
kVac = sqrt(kz.^2 - omega.^2);    
kCyl = sqrt(kz.^2 - epsilon * omega.^2);
       
% Determinant
I = besseli(n,kCyl);
Ip = besselip(n,kCyl);
K = besselk(n,kVac);
Kp = besselkp(n,kVac);
Delta = omega.^2 .* ((1./kVac).*(Kp./K) - epsilon.*(1./kCyl).*(Ip./I)) ...
    .* ((1./kVac).*(Kp./K) - (1./kCyl).*(Ip./I)) ...
    - n^2.*kz.^2.*(1./(kVac.^2) - 1./(kCyl.^2));

figure; hold on;
plot(kz, real(Delta), 'LineWidth', 1);
plot(kz, imag(Delta), 'LineWidth', 1);
ylim([-10, 10]);
xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\Delta/K_n$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting eigenmodes
close all;

kz = 2.372;
rhos = linspace(1e-3, 3, 120);
phi = linspace(-pi, pi, 240);

[R, P] = meshgrid(rhos, phi);
X = R.*cos(P); Y = R.*sin(P);

kVac = sqrt(kz.^2 - omega.^2);    
kCyl = sqrt(kz.^2 - epsilon * omega.^2);

I = besseli(n,kCyl);
Ip = besselip(n,kCyl);
K = besselk(n,kVac);
Kp = besselkp(n,kVac);

M11 = 1j*omega.*(besselkp(n,kVac)./kVac - epsilon.*besselip(n,kCyl)./besseli(n,kCyl).*besselk(n,kVac)./kCyl);
M12 = n.*kz.*(1./(kVac.^2) - 1./(kCyl.^2)) .* besselk(n,kVac);
M21 = M12;
M22 = -1j*omega.*(besselkp(n,kVac)./kVac - besselip(n,kCyl)./besseli(n,kCyl).*besselk(n,kVac)./kCyl);
    
Delta = M11 .* M22 - M12 .* M21;
disp(abs(Delta));

Bnk = 1;
eta0Dnk = - M11 / M12 * Bnk;

Ank = Bnk * K / I;
eta0Cnk = eta0Dnk * K / I;

EzFourier = zeros(1, numel(rhos));
EphiFourier = zeros(1, numel(rhos));
ErhoFourier = zeros(1, numel(rhos));
eta0HzFourier = zeros(1, numel(rhos));
eta0HphiFourier = zeros(1, numel(rhos));
eta0HrhoFourier = zeros(1, numel(rhos));

for i=1:numel(rhos)
    rho = rhos(i);
    EzFourier(i) = EzSecondaryFourier(Ank, Bnk, rho, n, kz, omega, epsilon);
    EphiFourier(i) = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    ErhoFourier(i) = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    eta0HzFourier(i) = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    eta0HphiFourier(i) = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
    eta0HrhoFourier(i) = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon);
end

EzFourier = real(exp(-1j*n*phi)'*EzFourier);
EphiFourier = real(exp(-1j*n*phi)'*EphiFourier);
ErhoFourier = real(exp(-1j*n*phi)'*ErhoFourier);
eta0HzFourier = real(exp(-1j*n*phi)'*eta0HzFourier);
eta0HphiFourier = real(exp(-1j*n*phi)'*eta0HphiFourier);
eta0HrhoFourier = real(exp(-1j*n*phi)'*eta0HrhoFourier);

I = EzFourier.^2 + 0*EphiFourier.^2 + 0*ErhoFourier.^2 + 0*eta0HzFourier.^2 ...
    + 0*eta0HphiFourier.^2 + 0*eta0HrhoFourier.^2;

figure; hold on;
surf(X, Y, I, 'EdgeColor', 'None');
plot3(cos(phi), sin(phi), 0.25 * ones(1, numel(phi)), '-w', 'LineWidth', 2);
view(2);
colorbar;

%% Plotting for all kz, omega
M = 1e3;
kz = linspace(-5, 5, M+1);
omega = linspace(0, 5, M+1);
[Kz, W] = meshgrid(kz, omega);

% Dispersion relation wavevectors
kVac = sqrt(Kz.^2 - W.^2);    
kCyl = sqrt(Kz.^2 - epsilon * W.^2);
       
% Determinant
I = besseli(n,kCyl);
Ip = besselip(n,kCyl);
K = besselk(n,kVac);
Kp = besselkp(n,kVac);
Delta = W.^2 .* ((1./kVac).*(Kp./K) - epsilon.*(1./kCyl).*(Ip./I)) ...
    .* ((1./kVac).*(Kp./K) - (1./kCyl).*(Ip./I)) ...
    - n^2.*Kz.^2.*(1./(kVac.^2) - 1./(kCyl.^2));

Delta = abs(Delta) < 0.25;

figure; hold on;
surf(Kz, W, real(Delta), 'EdgeColor', 'None');
% zlim([-5, 5]);
% caxis([-5, 5]);
xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
zlabel('$\Delta/K_n$', 'FontSize', 14, 'Interpreter', 'latex');
colorbar;