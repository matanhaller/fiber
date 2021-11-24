% Plotting field intensity in the time domain after calculation (iso-x
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
phi = linspace(-pi, pi, 240);
kz = linspace(-10, 10, 400);
% omega = [-2.822, 2.822];
omega = linspace(-10.0001, 10.0001, 400);
[K, W] = meshgrid(kz, omega);
N = -8:8;

dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);

y = [-flip(rhos), rhos];
z = linspace(-pi/dkz, pi/dkz, numel(kz)+1); z(end) = [];
z = z(151:250);
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
for t=215:225
    disp(t);
    
    EzInv = zeros(size(Y));
    EphiInv = zeros(size(Y));
    ErhoInv = zeros(size(Y));
    eta0HzInv = zeros(size(Y));
    eta0HphiInv = zeros(size(Y));
    eta0HrhoInv = zeros(size(Y));

    for i=1:numel(N)
        n = N(i);

        EzTotal = EzStructs(i).EzTotal(t-149,:,:);
        EzTotal = reshape(EzTotal, 100, numel(rhos));

        EphiTotal = EphiStructs(i).EphiTotal(t-149,:,:);
        EphiTotal = reshape(EphiTotal, 100, numel(rhos));

        ErhoTotal = ErhoStructs(i).ErhoTotal(t-149,:,:);
        ErhoTotal = reshape(ErhoTotal, 100, numel(rhos));

        eta0HzTotal = eta0HzStructs(i).eta0HzTotal(t-149,:,:);
        eta0HzTotal = reshape(eta0HzTotal, 100, numel(rhos));

        eta0HphiTotal = eta0HphiStructs(i).eta0HphiTotal(t-149,:,:);
        eta0HphiTotal = reshape(eta0HphiTotal, 100, numel(rhos));

        eta0HrhoTotal = eta0HrhoStructs(i).eta0HrhoTotal(t-149,:,:);
        eta0HrhoTotal = reshape(eta0HrhoTotal, 100, numel(rhos));

        EzInv(:,0.5*numel(y)+1:end) = EzInv(:,0.5*numel(y)+1:end) + (-1j)^n*EzTotal;
        EzInv(:,1:0.5*numel(y)) = EzInv(:,1:0.5*numel(y)) + (1j)^n*flip(EzTotal, 2);
        
        EphiInv(:,0.5*numel(y)+1:end) = EphiInv(:,0.5*numel(y)+1:end) + (-1j)^n*EphiTotal;
        EphiInv(:,1:0.5*numel(y)) = EphiInv(:,1:0.5*numel(y)) + (1j)^n*flip(EphiTotal, 2);
        
        ErhoInv(:,0.5*numel(y)+1:end) = ErhoInv(:,0.5*numel(y)+1:end) + (-1j)^n*ErhoTotal;
        ErhoInv(:,1:0.5*numel(y)) = ErhoInv(:,1:0.5*numel(y)) + (1j)^n*flip(ErhoTotal, 2);
        
        eta0HzInv(:,0.5*numel(y)+1:end) = eta0HzInv(:,0.5*numel(y)+1:end) + (-1j)^n*eta0HzTotal;
        eta0HzInv(:,1:0.5*numel(y)) = eta0HzInv(:,1:0.5*numel(y)) + (1j)^n*flip(eta0HzTotal, 2);
        
        eta0HphiInv(:,0.5*numel(y)+1:end) = eta0HphiInv(:,0.5*numel(y)+1:end) + (-1j)^n*eta0HphiTotal;
        eta0HphiInv(:,1:0.5*numel(y)) = eta0HphiInv(:,1:0.5*numel(y)) + (1j)^n*flip(eta0HphiTotal, 2);
        
        eta0HrhoInv(:,0.5*numel(y)+1:end) = eta0HrhoInv(:,0.5*numel(y)+1:end) + (-1j)^n*eta0HrhoTotal;
        eta0HrhoInv(:,1:0.5*numel(y)) = eta0HrhoInv(:,1:0.5*numel(y)) + (1j)^n*flip(eta0HrhoTotal, 2);
    end
    
    % Intensity (|E|^2 + |H|^2)
    I = real(EzInv).^2 + real(EphiInv).^2 + real(ErhoInv).^2 ...
      + real(eta0HzInv).^2 + real(eta0HphiInv).^2 + real(eta0HrhoInv).^2;

    % Poyting vector (E x H)
    Srho = real(EphiInv) .* real(eta0HzInv) - real(EzInv) .* real(eta0HphiInv);
    Sphi = - real(ErhoInv) .* real(eta0HzInv) + real(EzInv) .* real(eta0HrhoInv);
    Sz = real(ErhoInv) .* real(eta0HphiInv) - real(EphiInv) .* real(eta0HrhoInv);
    
    % Converting to Cartesian coordinates
%     Sy = Srho .* sign(Y);
    Sy = zeros(size(Y));
    
    figure; hold on;
    surf(Y, Z, log10(I), 'EdgeColor', 'none');
    quiver3(Y(1:5:end,1:5:end), Z(1:5:end,1:5:end), 10*ones(size(Y)/5), real(Sy(1:5:end,1:5:end)), real(Sz(1:5:end,1:5:end)), zeros(size(Y)/5), 'w', 'LineWidth', 1);
%     set(gcf, 'Visible', 'off');
    plot3(-ones(1, numel(z)), z, 10*ones(1, numel(z)), '-w', 'LineWidth', 2);
    plot3(ones(1, numel(z)), z, 10*ones(1, numel(z)), '-w', 'LineWidth', 2);
    view(2);
    xlim([-3,3]);
    ylim([-15,15]);
    xlabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$z$', 'FontSize', 14, 'Interpreter', 'latex');
    title(sprintf('Field Intensity: $\\varepsilon=%d$', epsilon), 'FontSize', 14, 'Interpreter', 'latex');
    caxis([-4,2]);
    colorbar;
%     saveas(gcf, sprintf('Simulation/tt=%d.jpg', t));
end
toc