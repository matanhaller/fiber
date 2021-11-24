% Plotting field intensity in the time domain after calculation (iso-z
% surfaces)

% Clearing workspace
close all;
clear;
clc;

epsilon = 4;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 0.8;

M = 120;
rhos = linspace(1e-3, 5, M);
phi = linspace(-pi, pi, 120);
kz = linspace(-10, 10, 400);
omega = linspace(-10.0001, 10.0001, 400);
[K, W] = meshgrid(kz, omega);
N = -8:8;

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);

[R, P] = meshgrid(rhos, phi);
X = R.*cos(P); Y = R.*sin(P);

% z = 0.5*numel(kz) + 1;
z = 51;

EzStructs = struct([]);
EphiStructs = struct([]);
ErhoStructs = struct([]);
eta0HzStructs = struct([]);
eta0HphiStructs = struct([]);
eta0HrhoStructs = struct([]);

%% Loading data
for i=1:numel(N)
    n = N(i);
    disp(n);
    if i == 1
        EzStructs = load(sprintf('Time Domain/epsilon=4,beta=0.8/Ez_n=%d.mat', n));
        EphiStructs = load(sprintf('Time Domain/epsilon=4,beta=0.8/Ephi_n=%d.mat', n));
        ErhoStructs = load(sprintf('Time Domain/epsilon=4,beta=0.8/Erho_n=%d.mat', n));
        eta0HzStructs = load(sprintf('Time Domain/epsilon=4,beta=0.8/eta0Hz_n=%d.mat', n));
        eta0HphiStructs = load(sprintf('Time Domain/epsilon=4,beta=0.8/eta0Hphi_n=%d.mat', n));
        eta0HrhoStructs = load(sprintf('Time Domain/epsilon=4,beta=0.8/eta0Hrho_n=%d.mat', n));
    else
        EzStructs(i) = load(sprintf('Time Domain/epsilon=4,beta=0.8/Ez_n=%d.mat', n));
        EphiStructs(i) = load(sprintf('Time Domain/epsilon=4,beta=0.8/Ephi_n=%d.mat', n));
        ErhoStructs(i) = load(sprintf('Time Domain/epsilon=4,beta=0.8/Erho_n=%d.mat', n));
        eta0HzStructs(i) = load(sprintf('Time Domain/epsilon=4,beta=0.8/eta0Hz_n=%d.mat', n));
        eta0HphiStructs(i) = load(sprintf('Time Domain/epsilon=4,beta=0.8/eta0Hphi_n=%d.mat', n));
        eta0HrhoStructs(i) = load(sprintf('Time Domain/epsilon=4,beta=0.8/eta0Hrho_n=%d.mat', n));
    end
end

%% Plotting intensity + Poynting vector

tic
for t=205:215
    disp(t);
    
    EzInv = zeros(size(R));
    EphiInv = zeros(size(R));
    ErhoInv = zeros(size(R));
    eta0HzInv = zeros(size(R));
    eta0HphiInv = zeros(size(R));
    eta0HrhoInv = zeros(size(R));

    for i=1:numel(N)
        n = N(i);
        
        EzTotal = EzStructs(i).EzTotal(t-149,z,:);
        EzTotal= EzTotal(:).';
        
        EphiTotal = EphiStructs(i).EphiTotal(t-149,z,:);
        EphiTotal = EphiTotal(:).';
        
        ErhoTotal = ErhoStructs(i).ErhoTotal(t-149,z,:);
        ErhoTotal = ErhoTotal(:).';

        eta0HzTotal = eta0HzStructs(i).eta0HzTotal(t-149,z,:);
        eta0HzTotal = eta0HzTotal(:).';
        
        eta0HphiTotal = eta0HphiStructs(i).eta0HphiTotal(t-149,z,:);
        eta0HphiTotal = eta0HphiTotal(:).';

        eta0HrhoTotal = eta0HrhoStructs(i).eta0HrhoTotal(t-149,z,:);
        eta0HrhoTotal = eta0HrhoTotal(:).';

        EzInv = EzInv + exp(-1j*n*phi).'*EzTotal;
        EphiInv = EphiInv + exp(-1j*n*phi).'*EphiTotal;
        ErhoInv = ErhoInv + exp(-1j*n*phi).'*ErhoTotal;
        eta0HzInv = eta0HzInv + exp(-1j*n*phi).'*eta0HzTotal;
        eta0HphiInv = eta0HphiInv + exp(-1j*n*phi).'*eta0HphiTotal;
        eta0HrhoInv = eta0HrhoInv + exp(-1j*n*phi).'*eta0HrhoTotal;
    end
    
    % Primary fields
    ExP = 0*ExPrimary(X, Y, 0, -pi/domega + t*(2*pi/(numel(omega)*domega)), x0, y0, z0, beta);
    EyP = 0*EyPrimary(X, Y, 0, -pi/domega + t*(2*pi/(numel(omega)*domega)), x0, y0, z0, beta);
    EzP = 0*EzPrimary(X, Y, 0, -pi/domega + t*(2*pi/(numel(omega)*domega)), x0, y0, z0, beta);
    eta0HyP = 0*eta0HyPrimary(X, Y, 0, -pi/domega + t*(2*pi/(numel(omega)*domega)), x0, y0, z0, beta);
    eta0HzP = 0*eta0HzPrimary(X, Y, 0, -pi/domega + t*(2*pi/(numel(omega)*domega)), x0, y0, z0, beta);
    
    ErhoP = ExP .* cos(P) + EyP .* sin(P);
    EphiP = -ExP .* sin(P) + EyP .* cos(P);
    eta0HrhoP = eta0HyP .* sin(P);
    eta0HphiP = eta0HyP .* cos(P);
    
    % Intensity (|E|^2 + |H|^2)
    I = real(EzInv + EzP).^2 + real(EphiInv + EphiP).^2 + real(ErhoInv + ErhoP).^2 ...
      + real(eta0HzInv + eta0HzP).^2 + real(eta0HphiInv + eta0HphiP).^2 + real(eta0HrhoInv + eta0HrhoP).^2;

    % Poyting vector (E x H)
    Srho = real(EphiInv) .* real(eta0HzInv) - real(EzInv) .* real(eta0HphiInv);
    Sphi = - real(ErhoInv) .* real(eta0HzInv) + real(EzInv) .* real(eta0HrhoInv);
    Sz = real(ErhoInv) .* real(eta0HphiInv) - real(EphiInv) .* real(eta0HrhoInv);
    
    % Converting to Cartesian coordinates
    Sx = Srho .* cos(P) - Sphi .* sin(P);
    Sy = Srho .* sin(P) + Sphi .* cos(P);
    
    % Angular momentum of em field
    Lz = X .* Sy - Y .* Sx;
    
    figure; hold on;
    plot3(cos(phi), sin(phi), 100 * ones(1, numel(phi)), '-w', 'LineWidth', 2);
    plot(x0 + beta*(-pi/domega + t*(2*pi/(numel(omega)*domega))), y0, 'wo', 'LineWidth', 2);
%     set(gcf, 'Visible', 'off');
    surf(X, Y, log10(I), 'EdgeColor', 'none');
    quiver3(X(1:5:end,1:5:end), Y(1:5:end,1:5:end), 10*ones(size(X)/5), real(Sx(1:5:end,1:5:end)), real(Sy(1:5:end,1:5:end)), zeros(size(X)/5), 'w', 'LineWidth', 1);
    view(2);
    xlim([-5,5]);
    ylim([-5,5]);
    xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('Angular Momentum: $\\varepsilon=%d$', epsilon), 'FontSize', 14, 'Interpreter', 'latex');
    caxis([-1,1]);
    colorbar;
%     saveas(gcf, sprintf('Simulation/t=%d.jpg', t));
end
toc