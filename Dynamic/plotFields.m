% Plotting EM fields due to the moving charge and cylinder

% Clearing workspace
close all;
clear;
clc;

set(gcf, 'Renderer', 'Painter');

x = linspace(-5, 5, 100);
y = linspace(-5, 5, 100);
[X, Y] = meshgrid(x, y);
z = 1;
t = 0;
x0 = 0; y0 = 2; z0 = 0;
betaVec = [0, 0.1, 0.5, 0.75, 0.95, 0.99];

figure(1);
tiledlayout(3, 2);
figure(2);
tiledlayout(3, 2);

% Plotting spatial fields
for beta=betaVec
    Ez = EzPrimary(X, Y, z, t, x0, y0, z0, beta);
    eta0Hz = eta0HzPrimary(X, Y, z, t, x0, y0, z0, beta);
    
    figure(1);
    nexttile; hold on;
    surf(X, Y, Ez, 'EdgeColor', 'none');
    view(2);
    title(sprintf('$\\beta = %.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
    
    figure(2);
    nexttile; hold on;
    surf(X, Y, eta0Hz, 'EdgeColor', 'none');
    view(2);
    title(sprintf('$\\beta = %.2f$', beta), 'FontSize', 14, 'Interpreter', 'latex');
end

% Plotting Fourier transform w.r.t angle
beta = 0.5;
N = 100;
rho = linspace(0, 5, N);
phi = linspace(-pi, pi, N);
[R, P] = ndgrid(rho, phi);
X = R.*cos(P); Y = R.*sin(P);
Ez = EzPrimary(X, Y, z, t, x0, y0, z0, beta);
EzFourier = fftshift(ifft(ifftshift(Ez,2),N,2),2);
n = linspace(-pi/(phi(2)-phi(1)), pi/(phi(2)-phi(1)), N+1); n(end) = [];

figure(3); tiledlayout(2, 2);
for ord=0:3
    nexttile; hold on;
    EzTrunc = 0.5 * (EzFourier(:,N/2+1+ord) .* exp(-ord*j*P) + EzFourier(:,N/2+1-ord) .* exp(ord*j*P));
    if n > 0
        EzTrunc = 2 * EzTrunc;
    end
    surf(X, Y, abs(EzTrunc), 'EdgeColor', 'None');
    view(2);
    colorbar;
end

eta0Hz = eta0HzPrimary(X, Y, z, t, x0, y0, z0, beta);
eta0HzFourier = fftshift(ifft(ifftshift(eta0Hz,2),N,2),2);

figure(4); tiledlayout(2, 2);
for ord=0:3
    nexttile; hold on;
    eta0HzTrunc = 0.5 * (eta0HzFourier(:,N/2+1+ord) .* exp(-ord*j*P) + eta0HzFourier(:,N/2+1-ord) .* exp(ord*j*P));
    if ord > 0
        eta0HzTrunc = 2 * eta0HzTrunc;
    end
    surf(X, Y, abs(eta0HzTrunc), 'EdgeColor', 'None');
    view(2);
    colorbar;
end

figure(5); hold on;
for ord=0:8
    plot(rho, abs(EzFourier(:,N/2+1+ord)), 'LineWidth', 1);
    xlabel('$\rho$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$|E_z(\rho,n)|$', 'FontSize', 14, 'Interpreter', 'latex');
end

figure(6); hold on;
for ord=0:8
    plot(rho, abs(eta0HzFourier(:,N/2+1+ord)), 'LineWidth', 1);
    xlabel('$\rho$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$|\eta_0 H_z(\rho,n)|$', 'FontSize', 14, 'Interpreter', 'latex');
end


% Calculating Fourier transform of electric field
beta = 0.8;
N = 1000;
kz = linspace(-20, 20, N+1); kz(end) = [];
omegas = linspace(0, 1, 9);
rho = 1;

figure(7); tiledlayout(3, 3);
for i=1:numel(omegas)
    omega = omegas(i);
    nexttile; hold on;
    for n=0:5
        EzFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
        plot(kz, abs(EzFourier), 'LineWidth', 1);
        xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|E_z(k_z,n)|$', 'FontSize', 14, 'Interpreter', 'latex');
    end
    title(sprintf('$\\bar{\\omega}=%.1f$', omega), 'FontSize', 14, 'Interpreter', 'latex');
end

figure(8); tiledlayout(3, 3);
omega = linspace(-10, 10, N+1); omega(end) = [];
for kz=0:8
    nexttile; hold on;
    for n=0:5
        EzFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
        plot(omega, abs(EzFourier), 'LineWidth', 1);
        xlabel('$\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|E_z(\bar{\omega},n)|$', 'FontSize', 14, 'Interpreter', 'latex');
    end
    title(sprintf('$k_z=%d$', kz), 'FontSize', 14, 'Interpreter', 'latex');
end

figure(9); tiledlayout(2, 2);
rho = linspace(-2, 2, N);
omega = 0;
for kz=0:3
    nexttile; hold on;
    for n=0:5
        EzFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
        plot(rho, abs(EzFourier), 'LineWidth', 1);
        xlabel('$\rho$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|E_z(\rho,n)|$', 'FontSize', 14, 'Interpreter', 'latex');
    end
    title(sprintf('$k_z=%d$', kz), 'FontSize', 14, 'Interpreter', 'latex');
end

% Calculating Fourier transform of magnetic field
N = 1000;
kz = linspace(-10, 10, N+1); kz(end) = [];
rho = 1;

figure(10); tiledlayout(2, 2);
for omega=0:3
    nexttile; hold on;
    for n=0:5
        eta0HzFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
        plot(kz, abs(eta0HzFourier), 'LineWidth', 1);
        xlabel('$k_z$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|\eta_0 H_z(k_z,n)|$', 'FontSize', 14, 'Interpreter', 'latex');
    end
    title(sprintf('$\\bar{\\omega}=%d$', omega), 'FontSize', 14, 'Interpreter', 'latex');
end

figure(11); tiledlayout(2, 2);
omega = linspace(-10, 10, N+1); omega(end) = [];
for kz=0:3
    nexttile; hold on;
    for n=0:5
        eta0HzFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
        plot(omega, abs(eta0HzFourier), 'LineWidth', 1);
        xlabel('$\bar{\omega}$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|\eta_0 H_z(\bar{\omega},n)|$', 'FontSize', 14, 'Interpreter', 'latex');
    end
    title(sprintf('$k_z=%d$', kz), 'FontSize', 14, 'Interpreter', 'latex');
end

figure(12); tiledlayout(2, 2);
rho = linspace(-2, 2, N);
omega = 0;
for kz=0:3
    nexttile; hold on;
    for n=0:5
        eta0HzFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
        plot(rho, abs(eta0HzFourier), 'LineWidth', 1);
        xlabel('$\rho$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$|\eta_0 H_z(\rho,n)|$', 'FontSize', 14, 'Interpreter', 'latex');
    end
    title(sprintf('$k_z=%d$', kz), 'FontSize', 14, 'Interpreter', 'latex');
end