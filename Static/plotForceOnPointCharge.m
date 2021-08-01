% Plotting the normalized force on the point charge

close all;
clear all;
clc;

% Plotting for different values of epsilonR
epsilonR = [1.2, 2, 4, 12];
eta = logspace(log10(1.02), log10(10), 100);
N = 80;
K = 600;

figure;
for er=epsilonR
    F = forceOnPointCharge(er, eta, N, 1e-2, K, 1e-3);
    F0 = -0.25 * (er - 1) / (er + 1) * 1 ./ ((eta - 1) .^ 2);
    Fnorm = F ./ F0;
    plot(eta - 1, Fnorm, 'LineWidth', 2, 'DisplayName', sprintf('$\\varepsilon=%.1f$', er)); hold on;
end
legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{F}$', 'FontSize', 14, 'Interpreter', 'latex');

% Plotting for different values of N
epsilonR = 2;
N = [0, 1, 2, 4, 10, 20, 40, 80];

figure; hold on;
for n=N
    F = forceOnPointCharge(epsilonR, eta, n, K, 1e-2, 1e-3);
    F0 = -0.25 * (er - 1) / (er + 1) * 1 ./ ((eta - 1) .^ 2);
    Fnorm = -F ./ F0;
    plot(eta - 1, Fnorm, 'LineWidth', 2, 'DisplayName', sprintf('$N=%d$', n));
end
legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{F}$', 'FontSize', 14, 'Interpreter', 'latex');

