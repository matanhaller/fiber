% Plotting the spectrum of radiation due to guided modes (in the far
% field regime).

close all;
clear;
clc;

N = -10:10;
epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 0;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = [0.5, 0.9, 1.1, 2] * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
nuMax = 40;

M = 1e3;
rhos = linspace(1e-3, 10, 120);
drho = rhos(2) - rhos(1);
kz = linspace(-10, 10, 40);
dkz = kz(2) - kz(1);
omegas = linspace(1e-2, 10.0001, 10);
[K, W] = meshgrid(kz, omegas);

%% Plotting for different velocities
tic;

figure;

for beta=betas
    spec = zeros(1, numel(omegas));
    
    for n=N
        disp(n);
        
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax);
        
        for rho=rhos
        
            EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, K, W, epsilon);
            
            EphiInv = fftshift(fft(fftshift(EphiFourier, 2), [], 2), 2);
            ErhoInv = fftshift(fft(fftshift(ErhoFourier, 2), [], 2), 2);
            eta0HphiInv = fftshift(fft(fftshift(eta0HphiFourier, 2), [], 2), 2);
            eta0HrhoInv = fftshift(fft(fftshift(eta0HrhoFourier, 2), [], 2), 2);
            
            EphiInv = EphiInv(:,end-1);
            ErhoInv = ErhoInv(:,end-1);
            eta0HphiInv = eta0HphiInv(:,end-1);
            eta0HrhoInv = eta0HrhoInv(:,end-1);
            
            Sz = conj(ErhoInv) .* eta0HphiInv - conj(EphiInv) .* eta0HrhoInv;
            spec = spec + rho * drho * Sz.';
        end
    end
    
    semilogy(omegas, real(spec), 'LineWidth', 1, 'DisplayName', sprintf('$\\beta=%.2f$', beta)); hold on;
end

toc;

xlabel('$\omega R/c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Normalized Spectrum', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');