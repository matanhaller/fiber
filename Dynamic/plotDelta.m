% Plotting determinant of the secondary fields' coefficient matrix.

% Clearing workspace
close all;
clear;
clc;

epsilon = 2; % Relative perimittivity
M = 1e3;
n = 0;

%% Plotting for kz=0 (Cutoff frequency)
kz = 0;
omega = linspace(-3, 3, M+1);

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
plot(omega, abs(Delta), 'LineWidth', 1);
ylim([0, 10]);
xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\Delta/K_n$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting for all kz, omega
M = 1e2;
kz = linspace(-3, 3, M+1);
omega = linspace(-3, 3, M+1);
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

figure; hold on;
surf(Kz, W, log10(abs(Delta)), 'EdgeColor', 'None');
xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
zlabel('$\Delta/K_n$', 'FontSize', 14, 'Interpreter', 'latex');
colorbar;