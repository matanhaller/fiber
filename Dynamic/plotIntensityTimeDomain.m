% Plotting field intensity in the time domain after calculation

% Clearing workspace
close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 0.8;

M = 120;
rhos = linspace(1e-3, 3, M);
phi = linspace(-pi, pi, 100);
kz = linspace(-8, 8, 40);
omega = linspace(-8.0001, 8.0001, 40);
[K, W] = meshgrid(kz, omega);

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);

[R, P] = meshgrid(rhos, phi);
X = R.*cos(P); Y = R.*sin(P);

z = 0.5*numel(kz) + 1;
for t=15:25
    EzInv = zeros(size(R));
    EphiInv = zeros(size(R));
    ErhoInv = zeros(size(R));
    eta0HzInv = zeros(size(R));
    eta0HphiInv = zeros(size(R));
    eta0HrhoInv = zeros(size(R));

    for n=-2:2

        EzTotal = load(sprintf('Time Domain/Ez_n=%d.mat', n));
        EzTotal = EzTotal.EzTotal;
        EzTotal = EzTotal(t,z,:);
        EzTotal= EzTotal(:)';

        EphiTotal = load(sprintf('Time Domain/Ephi_n=%d.mat', n));
        EphiTotal = EphiTotal.EphiTotal;
        EphiTotal = EphiTotal(t,z,:);
        EphiTotal = EphiTotal(:)';

        ErhoTotal = load(sprintf('Time Domain/Erho_n=%d.mat', n));
        ErhoTotal = ErhoTotal.ErhoTotal;
        ErhoTotal = ErhoTotal(t,z,:);
        ErhoTotal = ErhoTotal(:)';

        eta0HzTotal = load(sprintf('Time Domain/eta0Hz_n=%d.mat', n));
        eta0HzTotal = eta0HzTotal.eta0HzTotal;
        eta0HzTotal = eta0HzTotal(t,z,:);
        eta0HzTotal = eta0HzTotal(:)';

        eta0HphiTotal = load(sprintf('Time Domain/eta0Hphi_n=%d.mat', n));
        eta0HphiTotal = eta0HphiTotal.eta0HphiTotal;
        eta0HphiTotal = eta0HphiTotal(t,z,:);
        eta0HphiTotal = eta0HphiTotal(:)';

        eta0HrhoTotal = load(sprintf('Time Domain/eta0Hrho_n=%d.mat', n));
        eta0HrhoTotal = eta0HrhoTotal.eta0HrhoTotal;
        eta0HrhoTotal = eta0HrhoTotal(t,z,:);
        eta0HrhoTotal = eta0HrhoTotal(:)';

        EzInv = EzInv + exp(-1j*n*phi)'*EzTotal;
        EphiInv = EphiInv + exp(-1j*n*phi)'*EphiTotal;
        ErhoInv = ErhoInv + exp(-1j*n*phi)'*ErhoTotal;
        eta0HzInv = eta0HzInv + exp(-1j*n*phi)'*eta0HzTotal;
        eta0HphiInv = eta0HphiInv + exp(-1j*n*phi)'*eta0HphiTotal;
        eta0HrhoInv = eta0HrhoInv + exp(-1j*n*phi)'*eta0HrhoTotal;
    end

    I = abs(EzInv).^2 + abs(EphiInv).^2 + abs(ErhoInv).^2 ...
      + abs(eta0HzInv).^2 + abs(eta0HphiInv).^2 + abs(eta0HrhoInv).^2;

    figure; hold on;
    surf(X, Y, log10(I), 'EdgeColor', 'none');
    plot(cos(phi), sin(phi), '-w', 'LineWidth', 2);
    plot(x0 + beta*(-pi/domega + t*(2*pi/(numel(omega)*domega))), y0, 'wo', 'LineWidth', 2);
    view(2);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('Field Intensity: $\\varepsilon=%d$', epsilon), 'FontSize', 14, 'Interpreter', 'latex');
    caxis([-3,0]);
    colorbar;
end