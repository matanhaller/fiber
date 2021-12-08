% Plotting field intensity in the time domain after calculation (iso-z
% surfaces)

% Clearing workspace
close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 0;
beta = 0.99;

sigma = 0;
nuMax = 100;

M = 120;
rhos = linspace(1e-3, 5, M);
phi = linspace(-pi, pi, 4*M);
[R, P] = meshgrid(rhos, phi);

X = R .* cos(P);
Y = R .* sin(P);

Nk = 400;
Nw = 400;
kz = linspace(-10.0001, 10.0001, Nk);
omega = linspace(-20, 20, Nw);
[K, W] = meshgrid(kz, omega);

dkz = kz(2) - kz(1);
dw = omega(2) - omega(1);

z = linspace(-pi/dkz, pi/dkz, Nk+1); z(end) = [];
t = linspace(-pi/dw, pi/dw, Nw+1); t(end) = [];

N = -20:20;

EzStructs = struct([]);
EphiStructs = struct([]);
ErhoStructs = struct([]);
eta0HzStructs = struct([]);
eta0HphiStructs = struct([]);
eta0HrhoStructs = struct([]);

zi = floor(0.5 * Nk + 1);

%% Loading data
for i=1:numel(N)
    n = N(i);
    disp(n);
    if i == 1
        EzStructs = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Ez_n=%d.mat', n));
        EphiStructs = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Ephi_n=%d.mat', n));
        ErhoStructs = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Erho_n=%d.mat', n));
        eta0HzStructs = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hz_n=%d.mat', n));
        eta0HphiStructs = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hphi_n=%d.mat', n));
        eta0HrhoStructs = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hrho_n=%d.mat', n));
    else
        EzStructs(i) = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Ez_n=%d.mat', n));
        EphiStructs(i) = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Ephi_n=%d.mat', n));
        ErhoStructs(i) = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/Erho_n=%d.mat', n));
        eta0HzStructs(i) = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hz_n=%d.mat', n));
        eta0HphiStructs(i) = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hphi_n=%d.mat', n));
        eta0HrhoStructs(i) = load(sprintf('Time Domain/epsilon=2,beta=0.99/eta=1.4/eta0Hrho_n=%d.mat', n));
    end
end

%% Plotting intensity + Poynting vector
close all;

% Flags for plot configuration
saveFlag = 0;
quiverFlag = 0;
chargeDotFlag = 1;
primaryFieldFlag = 0;

tic
for ti=195:201
    disp(ti);
    
    EzInv = zeros(size(R));
    EphiInv = zeros(size(R));
    ErhoInv = zeros(size(R));
    eta0HzInv = zeros(size(R));
    eta0HphiInv  = zeros(size(R));
    eta0HrhoInv = zeros(size(R));

    for i=1:numel(N)
        n = N(i);
        
        EzTotal = EzStructs(i).EzTotal(:,ti);
        EzTotal= EzTotal(:).';
        
        EphiTotal = EphiStructs(i).EphiTotal(:,ti);
        EphiTotal = EphiTotal(:).';
        
        ErhoTotal = ErhoStructs(i).ErhoTotal(:,ti);
        ErhoTotal = ErhoTotal(:).';

        eta0HzTotal = eta0HzStructs(i).eta0HzTotal(:,ti);
        eta0HzTotal = eta0HzTotal(:).';
        
        eta0HphiTotal = eta0HphiStructs(i).eta0HphiTotal(:,ti);
        eta0HphiTotal = eta0HphiTotal(:).';

        eta0HrhoTotal = eta0HrhoStructs(i).eta0HrhoTotal(:,ti);
        eta0HrhoTotal = eta0HrhoTotal(:).';

        EzInv = EzInv + exp(-1j*n*phi).'*EzTotal;
        EphiInv = EphiInv + exp(-1j*n*phi).'*EphiTotal;
        ErhoInv = ErhoInv + exp(-1j*n*phi).'*ErhoTotal;
        eta0HzInv = eta0HzInv + exp(-1j*n*phi).'*eta0HzTotal;
        eta0HphiInv = eta0HphiInv + exp(-1j*n*phi).'*eta0HphiTotal;
        eta0HrhoInv = eta0HrhoInv + exp(-1j*n*phi).'*eta0HrhoTotal;
    end

    EzInv = real(EzInv);
    EphiInv = real(EphiInv);
    ErhoInv = real(ErhoInv);
    eta0HzInv = real(eta0HzInv);
    eta0HphiInv = real(eta0HphiInv);
    eta0HrhoInv = real(eta0HrhoInv);

    % Primary fields
    ExP = ExPrimary(X, Y, 0, t(ti), x0, y0, z0, beta);
    EyP = EyPrimary(X, Y, 0, t(ti), x0, y0, z0, beta);
    EzP = EzPrimary(X, Y, 0, t(ti), x0, y0, z0, beta);
    eta0HyP = eta0HyPrimary(X, Y, 0, t(ti), x0, y0, z0, beta);
    eta0HzP = eta0HzPrimary(X, Y, 0, t(ti), x0, y0, z0, beta);
    
    ErhoP = ExP .* cos(P) + EyP .* sin(P);
    EphiP = -ExP .* sin(P) + EyP .* cos(P);
    eta0HrhoP = eta0HyP .* sin(P);
    eta0HphiP = eta0HyP .* cos(P);
    
    % Poynting vector (E x H)
    Srho = EphiInv .* eta0HzInv - EzInv .* eta0HphiInv;
    Sphi = - ErhoInv .* eta0HzInv + EzInv .* eta0HrhoInv;
    Sz = ErhoInv .* eta0HphiInv - EphiInv .* eta0HrhoInv;

    % Adding primary components
    if primaryFieldFlag == 1
        EzInv = EzInv + EzP;
        EphiInv = EphiInv + EphiP;
        ErhoInv = ErhoInv + ErhoP;
        eta0HzInv = eta0HzInv + eta0HzP;
        eta0HphiInv = eta0HphiInv + eta0HphiP;
        eta0HrhoInv = eta0HrhoInv + eta0HrhoP;
    end

    % Intensity (|E|^2 + |H|^2)
%     I = eta0HzInv;
    I = EzInv.^2 + EphiInv.^2 + ErhoInv.^2 + eta0HzInv.^2 + eta0HphiInv.^2 + eta0HrhoInv.^2;  
    
    % Converting to Cartesian coordinates
    Sx = Srho .* cos(P) - Sphi .* sin(P);
    Sy = Srho .* sin(P) + Sphi .* cos(P);
    
    % Angular momentum of EM field
    Lz = X .* Sy - Y .* Sx;
    
    % Plotting animation
    figure; hold on;

    if saveFlag == 1
        set(gcf, 'Visible', 'off');
    end

    plot3(cos(phi), sin(phi), 100 * ones(1, numel(phi)), '--w', 'LineWidth', 1);
    
    if chargeDotFlag == 1
        plot3(x0 + beta*t(ti), y0, 100, 'wo', 'LineWidth', 2);
    end

    surf(X, Y, I, 'EdgeColor', 'none');

    if quiverFlag == 1
        quiver3(X(1:5:end,1:5:end), Y(1:5:end,1:5:end), 10*ones(size(X)/5), real(Sx(1:5:end,1:5:end)), real(Sy(1:5:end,1:5:end)), zeros(size(X)/5), 'w', 'LineWidth', 1);
    end

    view(2);
    xlim([-5,5]);
    ylim([-5,5]);
    xlabel('$x/R$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$y/R$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$ct=%.2f R$', t(ti)), 'FontSize', 14, 'Interpreter', 'latex');
    caxis([0,10]);
    colorbar;
%     colormap jet;

    if saveFlag == 1
        saveas(gcf, sprintf('Simulation/t=%d.png', ti));
    end
end
toc