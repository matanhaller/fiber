% Plotting coefficients of the secondary fields.

% Clearing workspace
close all;
clear;
clc;

N = 1000;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 1e-3;
epsilon = 2;

Nz = 100;
Ntheta = 101;

theta = linspace(-pi, pi, Ntheta);
rho = linspace(1e-2, 2.5, Nz);

omegas = linspace(1e-2, 0.5, 9);
kz = linspace(-4, 4, Nz);
dkz = kz(2) - kz(1);
n = linspace(-50, 50, Ntheta);

[kk, nn] = ndgrid(kz, n);

%% Plotting Fourier coefficients

kz = linspace(-2, 2, N);
omegas = linspace(1e-6, 1, 9);

figure(1); tiledlayout(3, 3);
figure(2); tiledlayout(3, 3);
for omega=omegas
    kVac = sqrt(kz.^2 - omega.^2);    
    kCyl = sqrt(kz.^2 - epsilon^2 * omega.^2);

    figure(1);
    nexttile; hold on;
    figure(2);
    nexttile; hold on;
    
    for n=0:4
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon);
        
        figure(1);
        plot(kz, abs(Ank .* besseli(n, kCyl*0.5)), 'LineWidth', 1);
        xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|A(n,k_z,\omega)|$', 'FontSize', 14, 'Interpreter', 'latex');
        title(sprintf('$\\bar{\\omega}=%.1f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
        
        figure(2);
        plot(kz, angle(Ank .* besseli(n, kCyl*0.5)), 'LineWidth', 1);
        xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$\angle A(n,k_z,\omega)$', 'FontSize', 14, 'Interpreter', 'latex');
        title(sprintf('$\\bar{\\omega}=%.1f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    end
end

kzs = 0:8; kzs(1) = 1e-6;
omega = linspace(-10, 10, N);

figure; tiledlayout(3, 3);

for kz=kzs
    nexttile; hold on;
    for n=0:4
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon);
        plot(omega, abs(Ank), 'LineWidth', 1);
        xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|A(n,k_z,\omega)|$', 'FontSize', 14, 'Interpreter', 'latex');
        title(sprintf('$k_z=%.1f$', kz), 'FontSize', 14, 'Interpreter', 'latex');
    end
end

%% Plotting electric field (frequency domain)

figure(1); tiledlayout(3, 3);
figure(2); tiledlayout(3, 3);
figure(3); tiledlayout(3, 3);
figure(4); tiledlayout(3, 3);
for omega=omegas
    kVac = sqrt(kk.^2 - omega.^2);    
    kCyl = sqrt(kk.^2 - epsilon^2 * omega.^2);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(nn, kk, omega, x0, y0, z0, beta, epsilon);

    ErhoInv = zeros(numel(rho), numel(theta));
    EphiInv = zeros(numel(rho), numel(theta));
    EzInv = zeros(numel(rho), numel(theta));
    
    for i=1:numel(rho)
        disp(i);
        r = rho(i);
        Erho = EzPrimaryFourier(r, nn, kk, omega, x0, y0, z0, beta) + ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        Ephi = EphiPrimaryFourier(r, nn, kk, omega, x0, y0, z0, beta) + EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        Ez = ErhoPrimaryFourier(r, nn, kk, omega, x0, y0, z0, beta) + EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        
        Erhoi = fftshift(fft2(fftshift(Erho))) * dkz;
        Ephii = fftshift(fft2(fftshift(Ephi))) * dkz;
        Ezi = fftshift(fft2(fftshift(Ez))) * dkz;
        
        ErhoInv(i,:) = Erhoi(Nz/2+1,:);
        EphiInv(i,:) = Ephii(Nz/2+1,:);
        EzInv(i,:) = Ezi(Nz/2+1,:);
    end
    
    I = abs(ErhoInv).^2 + abs(EphiInv).^2 + abs(EzInv).^2;
    
    [R, T] = ndgrid(rho, theta);
    X = R.*cos(T);
    Y = R.*sin(T);
    
    figure(1);
    nexttile; hold on;
    surf(X, Y, abs(ErhoInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(2);
    nexttile; hold on;
    surf(X, Y, abs(EphiInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(3);
    nexttile; hold on;
    surf(X, Y, abs(EzInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(4);
    nexttile; hold on;
    surf(X, Y, log10(abs(I)), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
end

%% Plotting magnetic field (frequency domain)

figure(5); tiledlayout(3, 3);
figure(6); tiledlayout(3, 3);
figure(7); tiledlayout(3, 3);
figure(8); tiledlayout(3, 3);
for omega=omegas
    kVac = sqrt(kk.^2 - omega.^2);    
    kCyl = sqrt(kk.^2 - epsilon^2 * omega.^2);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(nn, kk, omega, x0, y0, z0, beta, epsilon);

    eta0HrhoInv = zeros(numel(rho), numel(theta));
    eta0HphiInv = zeros(numel(rho), numel(theta));
    eta0HzInv = zeros(numel(rho), numel(theta));
    
    for i=1:numel(rho)
        r = rho(i);
        eta0Hrho = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        eta0Hphi = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        eta0Hz = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        
        eta0Hrhoi = fftshift(fft2(fftshift(eta0Hrho))) * dkz;
        eta0Hphii = fftshift(fft2(fftshift(eta0Hphi))) * dkz;
        eta0Hzi = fftshift(fft2(fftshift(eta0Hz))) * dkz;
        
        eta0HrhoInv(i,:) = eta0Hrhoi(Nz/2+2,:);
        eta0HphiInv(i,:) = eta0Hphii(Nz/2+2,:);
        eta0HzInv(i,:) = eta0Hzi(Nz/2+2,:);
    end
    
    I = abs(eta0HrhoInv).^2 + abs(eta0HphiInv).^2 + abs(eta0HzInv).^2;
    
    [R, T] = ndgrid(rho, theta);
    X = R.*cos(T);
    Y = R.*sin(T);
    
    figure(5);
    nexttile; hold on;
    surf(X, Y, abs(eta0HrhoInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(6);
    nexttile; hold on;
    surf(X, Y, abs(eta0HphiInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(7);
    nexttile; hold on;
    surf(X, Y, abs(eta0HzInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(8);
    nexttile; hold on;
    surf(X, Y, abs(I), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\omega=%.2f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
end

%% Plotting electric field (different velocities)

figure(9); tiledlayout(3, 3);
figure(10); tiledlayout(3, 3);
figure(11); tiledlayout(3, 3);
figure(12); tiledlayout(3, 3);

epsilon = 4;
omega = 2;
betas = [0.01, 0.1, 0.2, 0.4, 0.45, 0.55, 0.75, 0.9, 0.95];

for beta=betas
    kVac = sqrt(kk.^2 - omega.^2);    
    kCyl = sqrt(kk.^2 - epsilon^2 * omega.^2);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(nn, kk, omega, x0, y0, z0, beta, epsilon);

    ErhoInv = zeros(numel(rho), numel(theta));
    EphiInv = zeros(numel(rho), numel(theta));
    EzInv = zeros(numel(rho), numel(theta));
    
    for i=1:numel(rho)
        r = rho(i);
        Erho = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        Ephi = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        Ez = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        
        Erhoi = fftshift(fft2(fftshift(Erho)));
        Ephii = fftshift(fft2(fftshift(Ephi)));
        Ezi = fftshift(fft2(fftshift(Ez)));
        
        ErhoInv(i,:) = Erhoi(Nz/2+1,:);
        EphiInv(i,:) = Ephii(Nz/2+1,:);
        EzInv(i,:) = Ezi(Nz/2+1,:);
    end
    
    I = abs(ErhoInv).^2 + abs(EphiInv).^2 + abs(EzInv).^2;
    
    [R, T] = ndgrid(rho, theta);
    X = R.*cos(T);
    Y = R.*sin(T);
    
    figure(9);
    nexttile; hold on;
    surf(X, Y, abs(ErhoInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(10);
    nexttile; hold on;
    surf(X, Y, abs(EphiInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(11);
    nexttile; hold on;
    surf(X, Y, abs(EzInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(12);
    nexttile; hold on;
    surf(X, Y, log10(abs(I)), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
end

%% Plotting magnetic field (different velocities)

figure(13); tiledlayout(3, 3);
figure(14); tiledlayout(3, 3);
figure(15); tiledlayout(3, 3);
figure(16); tiledlayout(3, 3);

epsilon = 4;
omega = 2;
betas = [0.01, 0.1, 0.2, 0.4, 0.45, 0.55, 0.75, 0.9, 0.95];

for beta=betas
    kVac = sqrt(kk.^2 - omega.^2);    
    kCyl = sqrt(kk.^2 - epsilon^2 * omega.^2);
    
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(nn, kk, omega, x0, y0, z0, beta, epsilon);

    eta0HrhoInv = zeros(numel(rho), numel(theta));
    eta0HphiInv = zeros(numel(rho), numel(theta));
    eta0HzInv = zeros(numel(rho), numel(theta));
    
    for i=1:numel(rho)
        r = rho(i);
        eta0Hrho = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        eta0Hphi = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        eta0Hz = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, r, nn, kk, omega, x0, y0, z0, beta, epsilon);
        
        eta0Hrhoi = fftshift(ifft2(fftshift(eta0Hrho)));
        eta0Hphii = fftshift(ifft2(fftshift(eta0Hphi)));
        eta0Hzi = fftshift(ifft2(fftshift(eta0Hz)));
        
        eta0HrhoInv(i,:) = eta0Hrhoi(Nz/2+2,:);
        eta0HphiInv(i,:) = eta0Hphii(Nz/2+2,:);
        eta0HzInv(i,:) = eta0Hzi(Nz/2+2,:);
    end
    
    I = abs(eta0HrhoInv).^2 + abs(eta0HphiInv).^2 + abs(eta0HzInv).^2;
    
    [R, T] = ndgrid(rho, theta);
    X = R.*cos(T);
    Y = R.*sin(T);
    
    figure(13);
    nexttile; hold on;
    surf(X, Y, abs(eta0HrhoInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(14);
    nexttile; hold on;
    surf(X, Y, abs(eta0HphiInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(15);
    nexttile; hold on;
    surf(X, Y, abs(eta0HzInv), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
    
    figure(16);
    nexttile; hold on;
    surf(X, Y, abs(I), 'EdgeColor', 'none');
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\beta=%.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    colorbar;
end

