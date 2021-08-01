% Compare Fourier components of the primary field to the spatial-domain
% expressions

close all;
clear;
clc;

x0 = 0; y0 = 2; z0 = 0;
beta = 0.1;

rhos = linspace(1e-3, 1.5, 10);
phi = pi/2;
[R, P] = meshgrid(rhos, phi);
X = R.*cos(P); Y = R.*sin(P);
z = 0;
t = 0;

%% Plotting spatial-domain primary field intensity
ExSpatial = ExPrimary(X, Y, z, t, x0, y0, z0, beta);
EySpatial = EyPrimary(X, Y, z, t, x0, y0, z0, beta);
EzSpatial = EzPrimary(X, Y, z, t, x0, y0, z0, beta);
IESpatial = abs(ExSpatial).^2 + abs(EySpatial).^2 + abs(EzSpatial).^2;

eta0HxSpatial = 0;
eta0HySpatial = eta0HyPrimary(X, Y, z, t, x0, y0, z0, beta);
eta0HzSpatial = eta0HzPrimary(X, Y, z, t, x0, y0, z0, beta);
Ieta0HSpatial = abs(eta0HxSpatial).^2 + abs(eta0HySpatial).^2 + abs(eta0HzSpatial).^2;

ISpatial = IESpatial;

figure; hold on;
% plot(rhos, ISpatial, 'LineWidth', 1);

% figure; hold on;
% surf(X, Y, ISpatial, 'EdgeColor', 'none');
% view(2);

%% Plotting Fourier expansion of primary field intensity
rhos = linspace(1e-3, 1.9, 10);
phi = pi/2;
kz = linspace(-6.1, 6.1, 40);
omega = linspace(-4.1, 4.1, 40);

[K, W] = meshgrid(kz, omega);

EzInv = zeros(1, numel(rhos));
EphiInv = zeros(1, numel(rhos));
ErhoInv = zeros(1, numel(rhos));
eta0HzInv = zeros(1, numel(rhos));
eta0HphiInv = zeros(1, numel(rhos));
eta0HrhoInv = zeros(1, numel(rhos));

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);

for n=-1:1
    disp(n);
    EzTotal = zeros(1, numel(rhos));
    EphiTotal = zeros(1, numel(rhos));
    ErhoTotal = zeros(1, numel(rhos));
    eta0HzTotal = zeros(1, numel(rhos));
    eta0HphiTotal = zeros(1, numel(rhos));
    eta0HrhoTotal = zeros(1, numel(rhos));

    for i=1:numel(rhos)
        rho = rhos(i);
        EzFourier = EzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta);
        EphiFourier = EphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta);
        ErhoFourier = ErhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta);
        eta0HzFourier = eta0HzPrimaryFourier(rho, n, K, W, x0, y0, z0, beta);
        eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, K, W, x0, y0, z0, beta);
        eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, K, W, x0, y0, z0, beta);
        
        EzTotal(i) = EzTotal(i) + sum(sum(EzFourier)) * dkz * domega;
        EphiTotal(i) = EphiTotal(i) + sum(sum(EphiFourier)) * dkz * domega;
        ErhoTotal(i) = ErhoTotal(i) + sum(sum(ErhoFourier)) * dkz * domega;
        eta0HzTotal(i) = eta0HzTotal(i) + sum(sum(eta0HzFourier)) * dkz * domega;
        eta0HphiTotal(i) = eta0HphiTotal(i) + sum(sum(eta0HphiFourier)) * dkz * domega;
        eta0HrhoTotal(i) = eta0HrhoTotal(i) + sum(sum(eta0HrhoFourier)) * dkz * domega;
    
        disp(rho);
    end
    
    EzInv = EzInv + EzTotal * exp(1j*phi*n);
    EphiInv = EphiInv + EphiTotal * exp(1j*phi*n);
    ErhoInv = ErhoInv + ErhoTotal * exp(1j*phi*n);
    eta0HzInv = eta0HzInv + eta0HzTotal * exp(1j*phi*n);
    eta0HphiInv = eta0HphiInv + eta0HphiTotal * exp(1j*phi*n);
    eta0HrhoInv = eta0HrhoInv + eta0HrhoTotal * exp(1j*phi*n);
end

IEFourier = abs(EzInv).^2 + abs(EphiInv).^2 + abs(ErhoInv).^2;
Ieta0HFourier = abs(eta0HzInv).^2 + abs(eta0HphiInv).^2 + abs(eta0HrhoInv).^2;
IFourier = IEFourier + Ieta0HFourier;

% figure; hold on;
% surf(X, Y, IFourier, 'EdgeColor', 'none');
% view(2);

plot(rhos, IFourier, 'LineWidth', 1);

%% Plotting logarithmic relative error between intensities
% relErr = abs(ISpatial - IFourier) ./ ISpatial;
% 
% figure; hold on;
% surf(X, Y, relErr * 100, 'EdgeColor', 'none');
% colorbar;
% view(2)
