% Plotting the truncation error as function of N

close all;
clear;
clc;

% Plotting for different values of epsilonR
epsilonR = [1.2, 2, 4, 12];
Ninf = 80;
K = 500;

% Plotting for small eta
eta = linspace(1.02, 2, 100);
n = [0, 1, 2, 4, 10, 20, 30, 40, 50, 60, 70, 80];

figure;
for er=epsilonR
    Winf = potentialEnergyOnPointCharge(er, eta, Ninf, 1e-2, K, 1e-3);
    rmse = zeros(1, numel(n));
    for i=1:numel(n)
        W = potentialEnergyOnPointCharge(er, eta, n(i), 1e-2, K, 1e-3);
        rmse(i) = relRMSE(Winf, W);
    end
    semilogy(n, rmse, '--o', 'LineWidth', 2, 'DisplayName', sprintf('$\\varepsilon=%.1f$', er)); hold on;
end

legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$N$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('RMSE', 'FontSize', 14, 'Interpreter', 'latex');

% Plotting for large eta
eta = linspace(2, 100, 100);
n = [0, 1, 2, 4, 10, 20, 30, 40, 50, 60, 70, 80];

figure;
for er=epsilonR
    Winf = potentialEnergyOnPointCharge(er, eta, Ninf, 1e-2, K, 1e-3);
    rmse = zeros(1, numel(n));
    for i=1:numel(n)
        W = potentialEnergyOnPointCharge(er, eta, n(i), 1e-2, K, 1e-3);
        rmse(i) = relRMSE(Winf, W);
    end
    semilogy(n, rmse, '--o', 'LineWidth', 2, 'DisplayName', sprintf('$\\varepsilon=%.1f$', er)); hold on;
end

legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$N$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('RMSE', 'FontSize', 14, 'Interpreter', 'latex');
