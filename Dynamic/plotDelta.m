% Plotting determinant of the secondary fields' coefficient matrix.

% Clearing workspace
close all;
clear;
clc;

epsilon = 2; % Relative perimittivity
M = 1e4;
n = 1;

%% Plotting for kz=0 (Cutoff frequency)
close all;

kz = 3;
omega = linspace(0, 4.5, M+1);

% Dispersion relation wavevectors
kVac = sqrt(kz.^2 - omega.^2);    
kCyl = sqrt(kz.^2 - epsilon * omega.^2);
       
% Determinant
I = besseli(n,kCyl);
Ip = besselip(n,kCyl);
K = besselk(n,kVac);
Kp = besselkp(n,kVac);
Delta = DeltaCylinder(n, kz, omega, epsilon);

figure; hold on;
plot(omega, real(Delta), 'LineWidth', 1);
plot(omega, imag(Delta), 'LineWidth', 1);
plot(omega, zeros(1, numel(omega)), 'LineWidth', 1);
plot(kz/sqrt(epsilon)*ones(1,2), [-10,10], 'LineWidth', 1);
ylim([-10, 10]);
xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\Delta$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting for specific omega
close all;

n = 3;
omega = 8;
kz = linspace(0, 10, M+1);

Delta = DeltaCylinder(n, kz, omega, epsilon);

figure; hold on;
plot(kz, real(Delta), 'LineWidth', 1);
plot(kz, imag(Delta), 'LineWidth', 1);
plot(sqrt(epsilon)*omega*ones(1,2), [-0.1,0.1], 'LineWidth', 1);
ylim([-0.1, 0.1]);
xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\Delta$', 'FontSize', 14, 'Interpreter', 'latex');

% Plotting eigenmodes
kz = fzero(@(kz) real(DeltaCylinder(n, kz, omega, epsilon)), 8.5)
rhos = linspace(1e-3, 1.5, 120);
phi = linspace(-pi, pi, 1e3);

[R, P] = meshgrid(rhos, phi);
X = R.*cos(P); Y = R.*sin(P);

kVac = sqrt(kz.^2 - omega.^2); kVac = real(kVac) + 1j*sign(omega).*imag(kVac);    
kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);

I = besseli(n,kCyl);
Ip = besselip(n,kCyl);
K = besselk(n,kVac);
Kp = besselkp(n,kVac);
    
M11 = 1j*omega.*(epsilon.*Ip./kCyl - (Kp./K).*I./kVac) + 1e-100;
    M12 = n.*kz.*(1./(kCyl.^2) - 1./(kVac.^2)) .* I + 1e-100;
    M21 = M12;
    M22 = -1j*omega.*(Ip./kCyl - (Kp./K).*I./kVac) + 1e-100;

Delta = M11 .* M22 - M12 .* M21;

disp(abs(Delta));

Ank = 1;
Bnk = Ank * I / K;
eta0Cnk = - M11 / M12 * Ank;
eta0Dnk = eta0Cnk * I / K;

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

figure; hold on;
plot(rhos, abs(EzFourier).^2./max(abs(EzFourier).^2), 'LineWidth', 1);
xlabel('$\rho/R$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Amplitude', 'FontSize', 14, 'Interpreter', 'latex');

EzFourier = real(exp(-1j*n*phi).'*EzFourier);
EphiFourier = real(exp(-1j*n*phi).'*EphiFourier);
ErhoFourier = real(exp(-1j*n*phi).'*ErhoFourier);
eta0HzFourier = real(exp(-1j*n*phi).'*eta0HzFourier);
eta0HphiFourier = real(exp(-1j*n*phi).'*eta0HphiFourier);
eta0HrhoFourier = real(exp(-1j*n*phi).'*eta0HrhoFourier);

I = EzFourier.^2 + 0*EphiFourier.^2 + 0*ErhoFourier.^2 + 0*eta0HzFourier.^2 ...
    + 0*eta0HphiFourier.^2 + 0*eta0HrhoFourier.^2;

figure; hold on;
surf(X, Y, I./max(max(I)), 'EdgeColor', 'None');
plot3(cos(phi), sin(phi), 2 * ones(1, numel(phi)), '--w', 'LineWidth', 1);
view(2);
colorbar;

%% Plotting for all kz, omega
close all;

epsilon = 4;
n = -1;

M = 400;
kz = linspace(-10.0001, 10.0001, M);
omega = linspace(-20, 20, M);
[W, Kz] = meshgrid(omega, kz);

% Dispersion relation wavevectors
kVac = sqrt(Kz.^2 - W.^2); kVac = real(kVac) + 1j*sign(W).*imag(kVac);    
kCyl = sqrt(Kz.^2 - epsilon * W.^2); kCyl = real(kCyl) + 1j*sign(W).*imag(kCyl);
       
% Determinant
tic
Delta = DeltaCylinder(n, Kz, W, epsilon);
toc

% Delta = Delta .* (abs(W) > abs(Kz));

figure; hold on;
surf(W, Kz, log10(abs(Delta)), 'EdgeColor', 'none');
view(2);
xlim([-20,20]);
ylim([-10,10]);
% zlim([-5, 5]);
% caxis([-5, 5]);
xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
colorbar;