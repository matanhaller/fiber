% Plotting field intensity in the time domain after calculation (iso-x
% surfaces)

% Clearing workspace
close all;
clear;
clc;

epsilon = 2;
x0 = 0; y0 = 1.05; z0 = 0;
beta = 0.99;

sigma = 0;
nuMax = 100;

M = 120;
rhos = linspace(1e-3, 5, M);
phi = linspace(-pi, pi, M);
[R, P] = meshgrid(rhos, phi);

X = R .* cos(P);
Y = R .* sin(P);

Nk = 100;
Nw = 100;
kz = linspace(-10, 10, Nk);
omega = linspace(-20.0001, 20.0001, Nw);
[K, W] = meshgrid(kz, omega);

dkz = kz(2) - kz(1);
dw = omega(2) - omega(1);

z = linspace(-pi/dkz, pi/dkz, Nk+1); z(end) = [];
t = linspace(-pi/dw, pi/dw, Nw+1); t(1) = [];

N = -10:10;

y = [-flip(rhos), rhos];
Ny = numel(y);

[Y, Z] = meshgrid(y, z);

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
        EzStructs = load(sprintf('Time Domain/Ez_n=%d.mat', n));
        EphiStructs = load(sprintf('Time Domain/Ephi_n=%d.mat', n));
        ErhoStructs = load(sprintf('Time Domain/Erho_n=%d.mat', n));
        eta0HzStructs = load(sprintf('Time Domain/eta0Hz_n=%d.mat', n));
        eta0HphiStructs = load(sprintf('Time Domain/eta0Hphi_n=%d.mat', n));
        eta0HrhoStructs = load(sprintf('Time Domain/eta0Hrho_n=%d.mat', n));
    else
        EzStructs(i) = load(sprintf('Time Domain/Ez_n=%d.mat', n));
        EphiStructs(i) = load(sprintf('Time Domain/Ephi_n=%d.mat', n));
        ErhoStructs(i) = load(sprintf('Time Domain/Erho_n=%d.mat', n));
        eta0HzStructs(i) = load(sprintf('Time Domain/eta0Hz_n=%d.mat', n));
        eta0HphiStructs(i) = load(sprintf('Time Domain/eta0Hphi_n=%d.mat', n));
        eta0HrhoStructs(i) = load(sprintf('Time Domain/eta0Hrho_n=%d.mat', n));
    end
end

%% Plotting intensity + Poynting vector
tic
for ti=40:2:60
    disp(ti);
    
    EzInv = zeros(size(Y));
    EphiInv = zeros(size(Y));
    ErhoInv = zeros(size(Y));
    eta0HzInv = zeros(size(Y));
    eta0HphiInv = zeros(size(Y));
    eta0HrhoInv = zeros(size(Y));

    for i=1:numel(N)
        n = N(i);
        disp(n);
        EzTotal = EzStructs(i).EzTotal(ti,:,:);
        EzTotal = reshape(EzTotal, Nk, M);
        
        EphiTotal = EphiStructs(i).EphiTotal(ti,:,:);
        EphiTotal = reshape(EphiTotal, Nk, M);
        
        ErhoTotal = ErhoStructs(i).ErhoTotal(ti,:,:);
        ErhoTotal = reshape(ErhoTotal, Nk, M);

        eta0HzTotal = eta0HzStructs(i).eta0HzTotal(ti,:,:);
        eta0HzTotal = reshape(eta0HzTotal, Nk, M);

        eta0HphiTotal = eta0HphiStructs(i).eta0HphiTotal(ti,:,:);
        eta0HphiTotal = reshape(eta0HphiTotal, Nk, M);

        eta0HrhoTotal = eta0HrhoStructs(i).eta0HrhoTotal(ti,:,:);
        eta0HrhoTotal = reshape(eta0HrhoTotal, Nk, M);

        EzInv(:,0.5*Ny+1:end) = EzInv(:,0.5*Ny+1:end) + (-1j)^n*EzTotal;
        EzInv(:,1:0.5*Ny) = EzInv(:,1:0.5*Ny) + (1j)^n*flip(EzTotal, 2);
        
        EphiInv(:,0.5*Ny+1:end) = EphiInv(:,0.5*Ny+1:end) + (-1j)^n*EphiTotal;
        EphiInv(:,1:0.5*Ny) = EphiInv(:,1:0.5*Ny) + (1j)^n*flip(EphiTotal, 2);
        
        ErhoInv(:,0.5*Ny+1:end) = ErhoInv(:,0.5*Ny+1:end) + (-1j)^n*ErhoTotal;
        ErhoInv(:,1:0.5*Ny) = ErhoInv(:,1:0.5*Ny) + (1j)^n*flip(ErhoTotal, 2);
        
        eta0HzInv(:,0.5*Ny+1:end) = eta0HzInv(:,0.5*Ny+1:end) + (-1j)^n*eta0HzTotal;
        eta0HzInv(:,1:0.5*Ny) = eta0HzInv(:,1:0.5*Ny) + (1j)^n*flip(eta0HzTotal, 2);
       
        eta0HphiInv(:,0.5*Ny+1:end) = eta0HphiInv(:,0.5*Ny+1:end) + (-1j)^n*eta0HphiTotal;
        eta0HphiInv(:,1:0.5*Ny) = eta0HphiInv(:,1:0.5*Ny) + (1j)^n*flip(eta0HphiTotal, 2);
        
        eta0HrhoInv(:,0.5*Ny+1:end) = eta0HrhoInv(:,0.5*Ny+1:end) + (-1j)^n*eta0HrhoTotal;
        eta0HrhoInv(:,1:0.5*Ny) = eta0HrhoInv(:,1:0.5*Ny) + (1j)^n*flip(eta0HrhoTotal, 2);
    end

    EzInv = real(EzInv);
    EphiInv = real(EphiInv);
    ErhoInv = real(ErhoInv);
    eta0HzInv = real(eta0HzInv);
    eta0HphiInv = real(eta0HphiInv);
    eta0HrhoInv = real(eta0HrhoInv);
    
    % Intensity (|E|^2 + |H|^2)
    I = EzInv.^2 + EphiInv.^2 + ErhoInv.^2 + eta0HzInv.^2 + eta0HphiInv.^2 + eta0HrhoInv.^2;

    % Poyting vector (E x H)
    Srho = EphiInv .* eta0HzInv - EzInv .* eta0HphiInv;
    Sphi = - ErhoInv .* eta0HzInv + EzInv .* eta0HrhoInv;
    Sz = ErhoInv .* eta0HphiInv - EphiInv .* eta0HrhoInv;
    
    % Converting to Cartesian coordinates
%     Sy = Srho .* sign(Y);
    Sy = zeros(size(Y));
    
    figure; hold on;
    surf(Y, Z, I, 'EdgeColor', 'none');
%     quiver3(Y(1:5:end,1:5:end), Z(1:5:end,1:5:end), 10*ones(size(Y)/5), real(Sy(1:5:end,1:5:end)), real(Sz(1:5:end,1:5:end)), zeros(size(Y)/5), 'w', 'LineWidth', 1);
%     set(gcf, 'Visible', 'off');
    plot3(-ones(1, numel(z)), z, 10*ones(1, numel(z)), '--w', 'LineWidth', 1);
    plot3(ones(1, numel(z)), z, 10*ones(1, numel(z)), '--w', 'LineWidth', 1);
    view(2);
    xlim([-3,3]);
    ylim([-15,15]);
    xlabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$z$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('$\\bar{t}=%.2f$', t(ti)), 'FontSize', 14, 'Interpreter', 'latex');
    caxis([0,4]);
    colorbar;
%     saveas(gcf, sprintf('Simulation/tt=%d.jpg', ti));
end
toc