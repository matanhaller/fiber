% Plotting the spectrum of radiation due to guided modes (in the far
% field regime).

close all;
clear;
clc;

N = -10:10;
epsilon = 12;
x0 = 0; y0 = 1.4; z0 = 0;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = (0.1:0.1:2) * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
nuMax = 40;

M = 2.5e3;
maxModes = 15;
rhos = linspace(1e-3, 5, 120);
drho = rhos(2) - rhos(1);
omegas = repmat(linspace(1e-3, 10.0001, M), maxModes, 1);
W = zeros(1, numel(betas));
sigma = 0;

%% Plotting for different velocities
tic;

figure;

for k=1:numel(betas)
    beta = betas(k);
    
    spec = zeros(1, size(omegas,2));
    
    for n=N
        disp(n);
        
        ks = load(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, n));
        ks = ks.kz;
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
    
    W(k) = trapz(omegas(1,:), 0.5*(4*pi)^3*real(spec));
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

beta = 0.8;
    
currNnz = 0;
wsLoc = [];
   
for n=0:5
    disp(n);
        
    ks = load(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, n));
    ks = ks.kz;
        
    for m=1:size(ks,2)
        if nnz(ks(:,m)) > currNnz
            wsLoc = [wsLoc, m];
            currNnz = nnz(ks(:,m));
        end
    end
        
    spec = zeros(1, numel(wsLoc));
    ws = omegas(:,wsLoc);
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

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Calculating contribution of each harmonic
close all;

nuMax = 100;

epsilon = 12;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = [0.5, 0.9, 1.1, 2] * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));

W = zeros(numel(betas), numel(N));

tic;

figure;

for k=1:numel(betas)
    beta = betas(k);
      
    for n=N
        spec = zeros(1, size(omegas,2));
        
        disp(n);
        
        ks = load(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, n));
        ks = ks.kz;
        deltaDeriv = max(DeltaCylinderDeriv(n, ks, omegas, epsilon), 1e-100);
        
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffsTimesDelta(n, ks, omegas, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);
        
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
        
        W(k,n-N(1)+1) = trapz(omegas(1,:), 0.5*(4*pi)^3*real(spec));
    end
    
    semilogy(N, W(k,:), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.2f$', beta)); hold on;
end

toc;

xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Calculating contribution of each eigenmode
close all;

N = 0:10;
epsilon = 12;
beta = 0.8;
nuMax = 40;

x0 = 0; y0 = 1.4; z0 = 0;

M = 2.5e3;
maxModes = 15;
rhos = linspace(1e-3, 5, 2.5e3);
drho = rhos(2) - rhos(1);
omegas = linspace(1e-3, 10.0001, M);
W = zeros(numel(N), numel(maxModes));

tic;

for ni=N
    disp(ni);
    
    ks = load(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, ni));
    ks = ks.kz;
    
    for s=1:maxModes
        disp(s);
        disp(nnz(ks(s,:)));
        if nnz(ks(s,:)) == 0
            continue;
        end
        spec = zeros(1, M);
        
        for n=[-ni,ni]
            deltaDeriv = max(DeltaCylinderDeriv(n, ks(s,:), omegas, epsilon), 1e-100);
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffsTimesDelta(n, ks(s,:), omegas, x0, y0, z0, beta, epsilon, nuMax);
        
            for i=1:numel(rhos)
                rho = rhos(i);

                EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks(s,:), omegas, epsilon);
                ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks(s,:), omegas, epsilon);
                eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks(s,:), omegas, epsilon);
                eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, ks(s,:), omegas, epsilon);

                EphiModes = conj(EphiFourier)./deltaDeriv .* (ks(s,:) > 0);
                ErhoModes = conj(ErhoFourier)./deltaDeriv .* (ks(s,:) > 0);
                eta0HphiModes = eta0HphiFourier./deltaDeriv .* (ks(s,:) > 0);
                eta0HrhoModes = eta0HrhoFourier./deltaDeriv .* (ks(s,:) > 0);

                Sz = ErhoModes .* eta0HphiModes - EphiModes .* eta0HrhoModes;
                spec = spec + Sz * rho * drho;
            end
        end
        
        W(ni+1,s) = trapz(omegas, 0.5*(4*pi)^3*real(spec));
    end
    
    if ni == 0
        W(ni+1,:) = W(ni+1,:) / 2;
    end
    
end

figure; hold on;
for s=1:6
    plot(0:(size(W,1)-1), W(:,s), '--o', 'LineWidth', 1, 'DisplayName', sprintf('$s=%d$', s));
    set(gca, 'YScale', 'log');
    xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$W_{n,s}$', 'FontSize', 14, 'Interpreter', 'latex');
end
legend('FontSize', 14, 'Interpreter', 'latex');

toc;

%% Calculating radiated energy as function of distance from cylinder
close all;

N = 0:10;
epsilon = 2;
beta = 0.8;
nuMax = 160;

x0 = 0; z0 = 0;
y0Vec = logspace(log10(1.05), 1, 20);

M = 2.5e3;
maxModes = 15;
rhos = linspace(1e-3, 5, 120);
drho = rhos(2) - rhos(1);
omegas = linspace(1e-3, 10.0001, M);
W = zeros(1, numel(y0Vec));

tic;

%% Plotting for different distances from cylinder
y0Vec = logspace(log10(1.05), 1, 40);
epsilon = 12;
N = -10:10;
beta = 0.8;
nuMax = 40;
W = zeros(1, numel(y0Vec));

tic;

figure;

for k=1:numel(y0Vec)
    y0 = y0Vec(k);
    
    spec = zeros(1, size(omegas,2));
    
    for n=N
        disp(n);
        
        ks = load(sprintf('Dispersion Curve/eps=%d,n=%d.mat', epsilon, n));
        ks = ks.kz;
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
    
    W(k) = trapz(omegas(1,:), 0.5*(4*pi)^3*real(spec));
end

toc;

plot(y0Vec - 1, W, 'LineWidth', 1);
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');