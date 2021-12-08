% Plotting the normalized potential energy on the point charge

close all;
clear;
clc;

% Plotting for different values of epsilonR
epsilonR = [1.2, 2, 4, 12];
eta = logspace(log10(1.02), 1, 100);
N = 80;
K = 600;

% Winf = potentialEnergyOnPointChargeConducting(eta, 50, 0, K, 1e-3);

figure;
for er=epsilonR
    W0 = potentialEnergyDielectricHalfPlane(er, eta);
    W = potentialEnergyByRangeOfEta(er, eta, N, K, 1e-3);
    plot(eta - 1, W, 'LineWidth', 2, 'DisplayName', sprintf('$\\varepsilon=%.1f$', er)); hold on;
end

plot(eta - 1, Winf, 'LineWidth', 2, 'DisplayName', '$\varepsilon \rightarrow \infty$');

W0 = potentialEnergyDielectricHalfPlane(er, eta);
% plot(eta - 1, 0.5 * (eta-1).^(-1), '--', 'LineWidth', 2, 'DisplayName', 'Asymptotic');

legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');

% Plotting for different values of N
epsilonR = 2;
N = [0, 1, 2, 4, 10, 20, 40, 80];

figure; hold on;
for n=N
    W = potentialEnergyOnPointCharge(epsilonR, eta, n, 1e-2, K, 1e-3);
    W0 = potentialEnergyDielectricHalfPlane(epsilonR, eta);
    Wnorm = W ./ W0;
    plot(eta - 1, Wnorm, 'LineWidth', 2, 'DisplayName', sprintf('$N=%d$', n));
end
legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');

% Plotting coefficients
epsilonR = [1.2, 2, 4, 10, 100, 1000, 1e6, 1e9];
k = logspace(-15, 2, 100);
N = 0;
K = 600;

Winf = potentialEnergyOnPointChargeConductingCoeff(N, k, 1, 100);

figure;
for er=epsilonR
    W = potentialEnergyOnPointChargeCoeff(N, k, er, 100);
    loglog(k, W, 'LineWidth', 2, 'DisplayName', sprintf('$\\varepsilon=%.1f$', er)); hold on;
end

loglog(k, Winf, 'LineWidth', 2, 'DisplayName', '$\varepsilon \rightarrow \infty$');

legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$k$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');


