% Plotting the spectrum of radiation due to guided modes (in the far
% field regime).

close all;
clear;
clc;

N = -10:10;
epsilon = 12;
x0 = 0; y0 = 1.1; z0 = 0;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = [0.5, 0.9, 1.1, 2] * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
nuMax = 40;

M = 2.5e3;
maxModes = 15;
rhos = linspace(1e-3, 5, 120);
drho = rhos(2) - rhos(1);
omegas = repmat(linspace(1e-3, 10.0001, M), maxModes, 1);

%% Plotting for different velocities
tic;

figure;

for beta=betas
    spec = zeros(1, size(omegas,2));
    
    for n=-10:10
        disp(n);
        
        ks = load(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, n)).kz;
        deltaDeriv = max(DeltaCylinderDeriv(n, ks, omegas, epsilon), 1e-100);
        
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffsTimesDelta(n, ks, omegas, x0, y0, z0, beta, epsilon, nuMax);
        
        for i=1:numel(rhos)
            rho = rhos(i);
        
            EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, omegas, epsilon);
            ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, omegas, epsilon);
            eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, omegas, epsilon);
            eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, omegas, epsilon);
            
            EphiModes = sum(conj(EphiFourier)./deltaDeriv .* (ks > 0), 1);
            ErhoModes = sum(conj(ErhoFourier)./deltaDeriv .* (ks > 0), 1);
            eta0HphiModes = sum(eta0HphiFourier./deltaDeriv .* (ks > 0), 1);
            eta0HrhoModes = sum(eta0HrhoFourier./deltaDeriv .* (ks > 0), 1);
            
            Sz = ErhoModes .* eta0HphiModes - EphiModes .* eta0HrhoModes;
            spec = spec + Sz * rho * drho;
        end
    end
    
    semilogy(omegas(1,:), 0.5*(4*pi)^3*real(spec), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.2f$', beta)); hold on;
end

toc;

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Plotting for excitation frequencies
close all;

tic;

figure; hold on;

epsilon = 12;
y0 = 1.1;
rhos = linspace(1e-3, 5, 120);
beta = 0.96;
N = 0:10;
omegaFirst = [];

for n=0:10
    disp(n);
        
    ks = load(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, n));
    ks = ks.kz;
    
    currNnz = 0;
    wsLoc = [];
        
    for m=1:size(ks,2)
        if nnz(ks(:,m)) > currNnz
            wsLoc = [wsLoc, m];
            currNnz = nnz(ks(:,m));
        end
    end
        
    spec = zeros(1, numel(wsLoc));
    ws = omegas(:,wsLoc);
    if numel(ws) > 0
        omegaFirst = [omegaFirst, ws(1,1)];
    end
    ks = ks(:,wsLoc);
        
    deltaDeriv = max(DeltaCylinderDeriv(n, ks, ws, epsilon), 1e-100);
        
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffsTimesDelta(n, ks, ws, x0, y0, z0, beta, epsilon, nuMax);
        
    for i=1:numel(rhos)
        rho = rhos(i);
        
        EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, ws, epsilon);
        ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, ws, epsilon);
        eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, ws, epsilon);
        eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks, ws, epsilon);
            
        EphiModes = sum(conj(EphiFourier)./deltaDeriv .* (ks > 0), 1);
        ErhoModes = sum(conj(ErhoFourier)./deltaDeriv .* (ks > 0), 1);
        eta0HphiModes = sum(eta0HphiFourier./deltaDeriv .* (ks > 0), 1);
        eta0HrhoModes = sum(eta0HrhoFourier./deltaDeriv .* (ks > 0), 1);
            
        Sz = ErhoModes .* eta0HphiModes - EphiModes .* eta0HrhoModes;
        spec = spec + Sz * rho * drho;
    end
    
    stem(ws(1,:), 0.5*(4*pi)^3*real(spec), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
end
    
toc;

set(gca, 'YScale', 'log');
xlabel('$\omega_{n,s} R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

figure; hold on;
stem(0:(numel(omegaFirst)-1), omegaFirst, 'LineWidth', 1);
xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\omega_{n,1}$', 'FontSize', 14, 'Interpreter', 'latex');

%% Plotting for modes of given phase velocity
close all;

tic;

figure(1); hold on;
figure(2); hold on;

epsilon = 4;
y0 = 1.4;
rhos = linspace(1e-3, 5, 120);
drho = rhos(2) - rhos(1);
bC = 1/sqrt(epsilon);
gbC = bC / sqrt(1 - bC^2);
gb = 2*gbC;
beta = gb / sqrt(gb^2 + 1);
N = 0:10;
M = 2.5e3;
omega = linspace(1e-3, 50, M);
nuMax = 360;

for n=10
    disp(n);
        
    d = DeltaCylinder(n, omega/beta, omega, epsilon);
    sp = spmak(augknt(omega,2), real(d));
    w0 = fnzeros(sp);
    w0 = w0(1,:);
    k0 = w0 / beta;
    
    spec = zeros(1, numel(w0));
    
    deltaDeriv = max(DeltaCylinderDeriv(n, w0/beta, w0, epsilon), 1e-100);
    [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffsTimesDelta(n, w0/beta, w0, x0, y0, z0, beta, epsilon, nuMax);
        
    for i=1:numel(rhos)
        rho = rhos(i);
        
        EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, w0/beta, w0, epsilon);
        ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, w0/beta, w0, epsilon);
        eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, w0/beta, w0, epsilon);
        eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, w0/beta, w0, epsilon);
            
        EphiModes = conj(EphiFourier)./deltaDeriv;
        ErhoModes = conj(ErhoFourier)./deltaDeriv;
        eta0HphiModes = eta0HphiFourier./deltaDeriv;
        eta0HrhoModes = eta0HrhoFourier./deltaDeriv;
            
        Sz = ErhoModes .* eta0HphiModes - EphiModes .* eta0HrhoModes;
        spec = spec + Sz * rho * drho;
    end
    
    figure(1);
    stem(w0, 0.5*(4*pi)^3*real(spec), 'LineWidth', 1, 'DisplayName', sprintf('$n=%d$', n));
    
    w0lim = beta/sqrt(beta^2*epsilon-1) * (n*pi/2 + pi/4 - atan(1/epsilon ...
        * sqrt((beta^2*epsilon-1)/(1-beta^2))) + pi * (0:1:(ceil(numel(w0)/2)-1)));
    
    figure(2);
    stem(1:numel(w0), w0, 'o', 'LineWidth', 1, 'DisplayName', 'Numerical');
    stem(1:2:numel(w0), w0lim, 'o', 'LineWidth', 1, 'DisplayName', 'Kotanjyan');
end
    
toc;

figure(1);
set(gca, 'YScale', 'log');
xlabel('$\omega_{n,s} R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

figure(2);
xlabel('$s$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\omega_{n,s}$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');
