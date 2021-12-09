% Plotting the spectrum of radiation due to propagating waves (in the far
% field regime).

close all;
clear;
clc;

%% Plotting spectrum for different permittivities
% close all;

epsilonVec = [2, 4, 12];
x0 = 0; y0 = 1.4; z0 = 0;

M = 4e2;
omegaMax = 10;
omega = linspace(1e-3, omegaMax, M);
kz = linspace(1e-3, omegaMax, M);
[K, W] = meshgrid(kz, omega);
spec = zeros(numel(epsilonVec), numel(omega));
nuMax = 100;
sigma = 0;
eps = 1e-2;
tic

for j=1:numel(epsilonVec)
    epsilon = epsilonVec(j);
    disp(j);
  
    beta = 0.99;

    for n=-20:20
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0, eps);

        specInt = (abs(Bnk).^2 + abs(eta0Dnk).^2) ./ (W.^2 - K.^2 + eps) .* (abs(W) > abs(K));
        specInt(isnan(specInt)) = 0;
        spec(j,:) = spec(j,:) + omega .* trapz(kz, specInt, 2).';
    end
end

spec = 0.5 * (2*pi)^3 * real(spec);
% spec = medfilt1(spec, 3, [], 2);

figure; hold on;

for j=1:size(spec,1)
    plot(omega, spec(j,:), 'LineWidth', 1);
end

toc

xlabel('$\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\mathrm{d}\bar{W}/\mathrm{d}\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting spectrum for different velocities + contribution of each harmonic
% close all;

betaVec = [0.8, 0.9, 0.95, 0.97, 0.99];
x0 = 0; y0 = 1.05; z0 = 0;

N = -20:20;
epsilon = 2;

M = 4e2;
omegaMax = 20;
omega = linspace(1e-3, omegaMax, M);
kz = linspace(1e-3, omegaMax, M);
[K, W] = meshgrid(kz, omega);
spec = zeros(numel(betaVec), numel(omega));
WHarmonics = zeros(numel(betaVec), numel(N));
nuMax = 100;
sigma = 0;
eps = 1e-2;

tic

for j=1:numel(betaVec)
    beta = betaVec(j);
    disp(j);
  
    for k=1:numel(N)
        n = N(k);
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0, eps);

        specInt = (abs(Bnk).^2 + abs(eta0Dnk).^2) ./ (W.^2 - K.^2 + eps) .* (abs(W) > abs(K));
        specInt(isnan(specInt)) = 0;

        specHarmonic = omega .* trapz(kz, specInt, 2).';
        WHarmonics(j,k) = trapz(omega, specHarmonic);
        spec(j,:) = spec(j,:) + specHarmonic;
    end
end

toc

spec = 0.5 * (2*pi)^3 * real(spec);
WHarmonics = 0.5 * (2*pi)^3 * real(WHarmonics);

figure; hold on;

for j=1:size(spec,1)
    plot(omega, spec(j,:), 'LineWidth', 1);
end

xlabel('$\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\mathrm{d}\bar{W}/\mathrm{d}\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');

figure; hold on;

for j=1:size(WHarmonics,1)
    plot(N, WHarmonics(j,:), '--o', 'LineWidth', 1);
end

xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}_n$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting spectrum for different distances from the cylinder
% close all;

y0Vec = [1.05, 1.4, 2];
x0 = 0; z0 = 0;

N = -20:20;
epsilon = 4;
beta = 0.8;

M = 1e2;
omegaMax = 20;
omega = linspace(1e-3, omegaMax, M);
kz = linspace(1e-3, omegaMax, M);
[K, W] = meshgrid(kz, omega);
spec = zeros(numel(betaVec), numel(omega));
nuMax = 100;
sigma = 0;
eps = 1e-8;

tic

for j=1:numel(y0Vec)
    y0 = y0Vec(j);
    disp(j);
  
    for n=-20:20
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0, eps);

        specInt = (abs(Bnk).^2 + abs(eta0Dnk).^2) ./ (W.^2 - K.^2 + eps) .* (abs(W) > abs(K));
        specInt(isnan(specInt)) = 0;

        specHarmonic = omega .* trapz(kz, specInt, 2).';
        spec(j,:) = spec(j,:) + omega .* trapz(kz, specInt, 2).';
    end
end

toc

spec = 0.5 * (2*pi)^3 * real(spec);

figure; hold on;

for j=1:size(spec,1)
    plot(omega, spec(j,:), 'LineWidth', 1);
end

xlabel('$\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\mathrm{d}\bar{W}/\mathrm{d}\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting angular contribution
close all;

epsilon = 2;
betaVec = [0.2, 0.8, 0.99];
x0 = 0; y0 = 1.4; z0 = 0;
N = -20:20;
M = 1e2;
nuMax = 100;
rho = 100;
omegaMax = 20;
omega = linspace(1e-3, omegaMax, M);
kz = linspace(1e-3, omegaMax, M);
[K, W] = meshgrid(kz, omega);
Nphi = 1e3;
phi = linspace(-pi, pi, Nphi+1); phi(end) = [];
sigma = 0;

Wphi = zeros(numel(betaVec), Nphi);

eps = 0;

tic

for i=1:numel(betaVec)
    beta = betaVec(i);    
    EzFourier = zeros(numel(kz), numel(omega), numel(N));
    EphiFourier = zeros(numel(kz), numel(omega), numel(N));
    eta0HzFourier = zeros(numel(kz), numel(omega), numel(N));
    eta0HphiFourier = zeros(numel(kz), numel(omega), numel(N));
        
    for k=1:numel(N)
        n = N(k);
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0, eps);
        EzFourier(:,:,k) = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon, eps);
        EphiFourier(:,:,k) = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);
        eta0HzFourier(:,:,k) = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);
        eta0HphiFourier(:,:,k) = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon, eps);
    end
    
    EzFourier = fftshift(fft(EzFourier, Nphi, 3), 3);
    EphiFourier = fftshift(fft(EphiFourier, Nphi, 3), 3);
    eta0HzFourier = fftshift(fft(eta0HzFourier, Nphi, 3), 3);
    eta0HphiFourier = fftshift(fft(eta0HphiFourier, Nphi, 3), 3);
        
    S = conj(EphiFourier) .* eta0HzFourier - conj(EzFourier) .* eta0HphiFourier;
    for k=1:numel(N)
        S(:,:,k) = S(:,:,k) .* (abs(W) > abs(K));
    end
    S(isnan(S)) = 0;
    
    Wphi(i,:) = 4*pi*rho * trapz(omega, trapz(kz, S, 2), 1);
end

toc

disp(max(abs(real(Wphi)), [], 2));
disp(max(abs(imag(Wphi)), [], 2));

Wphi = real(Wphi);

figure; hold on;

for i=1:numel(betaVec)
    plot(phi * (180/pi), Wphi(i,:), 'LineWidth', 1); 
end

xlim([-180, 180]);
xlabel('$\phi$ [$^{\circ}$]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\mathrm{d}\bar{W}/\mathrm{d}\phi$', 'FontSize', 14, 'Interpreter', 'latex');

figure;
polarplot(phi, Wphi ./ max(Wphi,[],2), 'LineWidth', 2);


%% Plotting angle-frequency spectrum
close all;

epsilon = 2;
beta = 0.99;
x0 = 0; y0 = 1.05; z0 = 0;
N = -20:20;
M = 100;
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
    
    for k=1:numel(N)
        n = N(k);
        disp(n);
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);
        EzFourier(:,:,k) = EzSecondaryFourier(Ank, Bnk, rho, n, K, W, epsilon);
        EphiFourier(:,:,k) = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0HzFourier(:,:,k) = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
        eta0HphiFourier(:,:,k) = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
    end

    EzFourier = fftshift(fft(EzFourier, Nphi, 3), 3);
    EphiFourier = fftshift(fft(EphiFourier, Nphi, 3), 3);
    eta0HzFourier = fftshift(fft(eta0HzFourier, Nphi, 3), 3);
    eta0HphiFourier = fftshift(fft(eta0HphiFourier, Nphi, 3), 3);
    
    S = conj(EphiFourier) .* eta0HzFourier - conj(EzFourier) .* eta0HphiFourier;
    for k=1:numel(N)
        S(:,:,k) = tril(S(:,:,k));
    end
    S(isnan(S)) = 0;

    Wphi = 4*pi*rho * trapz(kz, S, 2);
    Wphi = reshape(Wphi, numel(omega), Nphi);

%     disp(max(abs(real(Wphi))));
%     disp(max(abs(imag(Wphi))));

    surf(phi * (180/pi), omega, real(Wphi), 'EdgeColor', 'none');
    view(2);
end
toc

xlim([-180, 180]);
ylim([0, 12]);
xlabel('$\phi$ [$^{\circ}$]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
