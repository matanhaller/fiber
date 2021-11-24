% Plotting the spectrum of radiation due to propagating waves (in the far
% field regime).

close all;
clear;
clc;

N = -10:10;
epsilon = 12;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 0.49;
nuMax = 40;

M = 2.5e3;
rho = 1.5;
kz = linspace(1e-4, 0.75, M);
omegas = linspace(1e-2, 6.0001, M);
W = zeros(numel(epsilon), numel(omegas));

%% Plotting spectrum for different permittivities
tic
figure;
for j=1:numel(epsilon)
    er = epsilon(j);
    disp(er);
    for n=N
        disp(n);
        for i=1:numel(omegas)
            omega = omegas(i);
            disp(omega);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [K, Wi] = meshgrid(kRange, omega);
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, Wi, x0, y0, z0, beta, er, nuMax);
            W(j, i) = W(j, i) + omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        end
    end
end

for j=1:numel(epsilon)
    er = epsilon(j);
    semilogy(omegas, W(j, :), 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%d$', er)); hold on;
end

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');
toc

%% Plotting spectrum for different velocities
epsilon = 4;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = [0.5, 0.9, 1.1, 2] * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
W = zeros(numel(betas), numel(omegas));

tic
figure; hold on;
for j=1:numel(betas)
    beta = betas(j);
    disp(beta);
    for n=N
        disp(n);
        for i=1:numel(omegas)
            omega = omegas(i);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kRange, omega, x0, y0, z0, beta, epsilon);
            W(j, i) = W(j, i) + omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        end
    end
end

for j=1:numel(betas)
    beta = betas(j);
    plot(omegas, W(j, :) / max(W(j, :)), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.2f$', beta));
end

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');
toc

%% Integrating over spectrum for different harmonics
epsilon = 12;
N = -10:10;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = [0.5, 1, 2, 4] * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
Energy = zeros(numel(betas), numel(N));

tic
figure;
for j=1:numel(betas)
    beta = betas(j);
    disp(beta);
    for k=1:numel(N)
        n = N(k);
        disp(n);
        W = zeros(1, numel(omegas));
        for i=1:numel(omegas)
            omega = omegas(i);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kRange, omega, x0, y0, z0, beta, epsilon);
            W(i) = omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        end
        Energy(j, k) = trapz(omegas, W);
    end
end

for j=1:numel(betas)
    beta = betas(j);
    semilogy(N, Energy(j, :) ./ max(Energy(j, :)), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.2f$', beta)); hold on;
end

toc

xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Energy', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Plotting spectrum for different impact parameters
y0Vec = 10;
% y0Vec = [1.05, 1.4, 2, 5, 10];
epsilon = 12;

tic
figure;
for j=1:numel(y0Vec)
    y0 = y0Vec(j);
    for n=N
        disp(n);
        for i=1:numel(omegas)
            omega = omegas(i);
            disp(omega);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [K, Wi] = meshgrid(kRange, omega);
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, Wi, x0, y0, z0, beta, epsilon, nuMax);
            W(j, i) = W(j, i) + omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        end
    end
end

for j=1:numel(y0Vec)
    y0 = y0Vec(j);
    semilogy(omegas, W(j, :), 'LineWidth', 1, 'DisplayName', sprintf('$y_0=%d$', y0)); hold on;
end

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');
toc

%% Plotting spectrum for different distances from the cylinder
epsilon = 12;
beta = 0.8;
y0Vec = logspace(log10(1.05), 1, 40);
M = 400;
omega = linspace(0.0101, 12.0001, M);
domega = omega(2) - omega(1);
kz = linspace(1e-2, 12, M);
[K, Wi] = meshgrid(kz, omega);
W = zeros(1, numel(y0Vec));
spec = zeros(numel(y0Vec), numel(omega));

tic
figure;
for j=1:numel(y0Vec)
    y0 = y0Vec(j);
    disp(j);
    for n=N
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, Wi, x0, y0, z0, beta, epsilon, nuMax);
        specInt = tril(1./(Wi.^2 - K.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        spec(j,:) = spec(j,:) + omega .* trapz(domega, specInt, 2).';
    end
    W(j) = 0.5*(2*pi)^3 * trapz(omega, spec(j,:));
end
toc

plot(y0Vec - 1, W, 'LineWidth', 1);
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');
savefig(gcf, 'multiple_eta_eps_12.fig');

figure; hold on;
for i=1:numel(y0Vec)
    plot(omega, spec(i,:), 'LineWidth', 1);
end

savefig(gcf, 'spec_multiple_eta_eps_12.fig');

%% Plotting spectrum for different velocities (2nd attempt)
epsilon = 12;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = (0.1:0.1:4) * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
y0 = 1.4;

M = 400;
omega = linspace(0.0101, 12.0001, M);
domega = omega(2) - omega(1);
kz = linspace(1e-2, 12, M);
[K, Wi] = meshgrid(kz, omega);
W = zeros(1, numel(betas));
spec = zeros(numel(betas), numel(omega));

tic
figure;
for j=1:numel(betas)
    beta = betas(j);
    disp(j);
    for n=N
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, Wi, x0, y0, z0, beta, epsilon, nuMax);
        specInt = tril(1./(Wi.^2 - K.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
        spec(j,:) = spec(j,:) + omega .* trapz(domega, specInt, 2).';
    end
    W(j) = 0.5*(2*pi)^3 * trapz(omega, spec(j,:));
end
toc

plot(gb, W, 'LineWidth', 1);
xlabel('$\gamma\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');
savefig(gcf, 'multiple_beta_eps_12.fig');


figure; hold on;
for i=1:numel(betas)
    plot(omega, spec(i,:), 'LineWidth', 1);
end

savefig(gcf, 'spec_multiple_beta_eps_12.fig');

%% Plotting angular contribution
close all;

epsilon = 4;
beta = 0.8;
x0 = 0; y0 = 1.4; z0 = 0;
N = -10:10;
M = 400;
nuMax = 100;
rhos = 40;
omega = linspace(0.0101, 12.0001, M);
kz = linspace(1e-2, 12, M);
[K, W] = meshgrid(kz, omega);
Nphi = 1e3;
phi = linspace(-pi, pi, Nphi+1); phi(end) = [];
sigma = 0;

tic
figure; hold on;
for j=1:numel(rhos)
    rho = rhos(j);
%     disp(rho);

    EzFourier = zeros(numel(kz), numel(omega), numel(N));
    EphiFourier = zeros(numel(kz), numel(omega), numel(N));
    eta0HzFourier = zeros(numel(kz), numel(omega), numel(N));
    eta0HphiFourier = zeros(numel(kz), numel(omega), numel(N));

    Wphi = zeros(1, Nphi);
    
    for k=1:numel(N)
        n = N(k);
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);
        EzFourier(:,:,k) = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
        EphiFourier(:,:,k) = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0HzFourier(:,:,k) = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0HphiFourier(:,:,k) = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
    end

    EzFourier = fftshift(fft(EzFourier, Nphi, 3));
    EphiFourier = fftshift(fft(EphiFourier, Nphi, 3));
    eta0HzFourier = fftshift(fft(eta0HzFourier, Nphi, 3));
    eta0HphiFourier = fftshift(fft(eta0HphiFourier, Nphi, 3));
    
    S = conj(EphiFourier) .* eta0HzFourier - conj(EzFourier) .* eta0HphiFourier;
    for k=1:numel(N)
        S(:,:,k) = tril(S(:,:,k));
    end

    Wphi = 4*pi*rho * trapz(omega, trapz(kz, S, 2), 1);
    Wphi = Wphi(:);

    disp(max(abs(real(Wphi))));
    disp(max(abs(imag(Wphi))));

    plot(phi * (180/pi), real(Wphi), 'LineWidth', 1);
end
toc

xlabel('$\phi$ [$^{\circ}$]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\mathrm{d}\bar{W}/\mathrm{d}\phi$', 'FontSize', 14, 'Interpreter', 'latex');