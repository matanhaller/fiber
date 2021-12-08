% Plotting the work done by the cylinder's field on the electron as a
% function of its velocity, distance from the cylinder and the cylinder's
% permittivity.

close all;
clear;
clc;

N = -10:10;
epsilon = 2;
x0 = 0; y0 = 1.4; z0 = 0;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = (0.1:0.1:2) * gbC; 
betas = sqrt(gb.^2 ./ (gb.^2 + 1));
nuMax = 100;
sigma = 0;

xMax = 10;
dx = 1e-2;
x = -xMax:dx:xMax;
rho = sqrt(x.^2 + y0^2);
phi = atan(y0./x);

M = 40;
kz = linspace(-10, 10, M);
omega = linspace(-10.0001, 10.0001, M);
dkz = kz(2) - kz(1);
domega = omega(2) - omega(1);
[K, W] = meshgrid(kz, omega);

%% Plotting for different velocities
tic;

work = zeros(1, numel(betas));

for i=1:numel(betas)
    beta = betas(i);
    disp(beta);
    t = x / beta;
    
    for n=N
        disp(n);
        
        [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, K, W, x0, y0, z0, beta, epsilon, nuMax, sigma, 0);
        Ex = zeros(1, numel(x));
        
        for m=1:numel(x)
   
            EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho(m), n, K, W, epsilon);
            ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho(m), n, K, W, epsilon);
            
            ExFourier = ErhoFourier * x(m)/rho(m) - EphiFourier * y0/rho(m);
            
            ExInv = sum(ExFourier, 2) * dkz;
            ExInv = ExInv.';
            Ex(m) = trapz(omega, ExInv .* exp(1j*omega*t(m)));
              
        end
        
        work(i) = work(i) - trapz(x, Ex .* exp(-1j*n*phi));
    end
    
end

toc;

figure; hold on;

plot(gb / gbC, real(work), 'LineWidth', 1);
xlabel('$\gamma\beta/(\gamma\beta)_\mathrm{C}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Radiated Energy', 'FontSize', 14, 'Interpreter', 'latex');