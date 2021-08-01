% Asymptotic analysis of the potential energy on the point charge due to
% the fiber

close all;
clear;
clc;

% Plotting for different values of epsilonR
epsilonR = 12;
eta = logspace(log10(1.02), log10(1e3), 1e3);
N = 80;
K = 500;

figure;

W = potentialEnergyByRangeOfEta(epsilonR, eta, N, K, 1e-3);
loglog(eta - 1, W, 'LineWidth', 2, 'DisplayName', 'Exact');
hold on;

% for n=0:1
%     Wtrunc = potentialEnergyByRangeOfEta(epsilonR, eta, n, K, 1e-3);
%     loglog(eta - 1, Wtrunc, 'LineWidth', 2, 'DisplayName', sprintf('$N=%d$', n));
% end

% Dielectric half-plane potential energy (eta-->1)
W0 = potentialEnergyDielectricHalfPlane(epsilonR, eta);
loglog(eta - 1, W0, '--', 'LineWidth', 2, 'DisplayName', 'Asymptotic $\propto (\eta -1)^{-1}$');

% Induced dipole with ~1/(eta-1) dipole moment potential energy (eta-->inf)
a = 0.589;
b = 0.098;
Winf = ((epsilonR - 1) / (epsilonR + 1) * a + (epsilonR - 1) * b) * (eta - 1).^(-3);
loglog(eta - 1, Winf, '--', 'LineWidth', 2, 'DisplayName', 'Asymptotic $\propto  (\eta -1)^{-3}$');

p = 0.5;
Wtrans = (W0.^(-p) + Winf.^(-p)).^(-1/p);
loglog(eta - 1, Wtrans, '--', 'LineWidth', 2, 'DisplayName', 'Transient');

xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$W$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

% Checking various "transient" approximations
figure;

pVec = 0.4:1e-2:0.6;

rmseSmallEta = zeros(1, numel(pVec));
rmseLargeEta = zeros(1, numel(pVec));

for i=1:numel(pVec)
    p = pVec(i);
    Wtrans = (W0.^(-p) + Winf.^(-p)).^(-1/p);
    rmseSmallEta(i) = relRMSE(W(eta <= 10), Wtrans(eta <= 10));
    rmseLargeEta(i) = relRMSE(W(eta > 10), Wtrans(eta > 10));
    loglog(eta - 1, relError(W, Wtrans), 'LineWidth', 2, 'DisplayName', sprintf('$p=%.2f$', p));
    hold on;
end

loglog(eta - 1, relError(W, W0), 'LineWidth', 2, 'DisplayName', 'Asymptotic $\propto  (\eta -1)^{-1}$');
loglog(eta - 1, relError(W, Winf), 'LineWidth', 2, 'DisplayName', 'Asymptotic $\propto  (\eta -1)^{-3}$');

xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Error', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

% Plotting RMSE as a function of exponent
figure; hold on;

plot(pVec, rmseSmallEta, '--o', 'LineWidth', 2, 'DisplayName', '$\eta<10$');
plot(pVec, rmseLargeEta, '--o', 'LineWidth', 2, 'DisplayName', '$\eta>10$');

xlabel('$p$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('RMSE', 'FontSize', 14, 'Interpreter', 'latex');

legend('FontSize', 14, 'Interpreter', 'latex');

% Plotting relative error for various values of epsilon
% epsilonR = [1.2, 2, 4, 12];
% eta = logspace(log10(1.02), log10(1e3), 100);
% 
% figure;
% for er=epsilonR
%     W = potentialEnergyByRangeOfEta(er, eta, N, K, 1e-3);
%     Winf = ((er - 1) / (er + 1) * a + (er - 1) * b) * (eta - 1).^(-3);
%     
%     loglog(eta - 1, relError(W, Winf), 'LineWidth', 2, 'DisplayName', sprintf('$\\varepsilon=%.1f$', er)); hold on;
%     xlabel('$\eta - 1$', 'FontSize', 14, 'Interpreter', 'latex');
%     ylabel('Error', 'FontSize', 14, 'Interpreter', 'latex');
%     legend('FontSize', 14, 'Interpreter', 'latex');
% end

% Plotting for variable permittivity
epsilon = logspace(log10(1.01), log10(12), 40);
eta = 2.5;
Weps = zeros(1, numel(epsilon));

for i=1:numel(epsilon)
    er = epsilon(i);
    Weps(i) = potentialEnergyByRangeOfEta(er, eta, N, K, 1e-3);
end

[fitresult, gof] = createFit(epsilon, Weps);
disp(gof);

a = fitresult.a;
b = fitresult.b;

Wfit = a * (epsilon - 1) + b * (epsilon - 1) ./ (epsilon + 1);

figure; hold on;
plot(epsilon, Weps, 'LineWidth', 2, 'DisplayName', 'Numeric');
plot(epsilon, Wfit, '--', 'LineWidth', 2, 'DisplayName', 'Analytic');
xlabel('$\varepsilon$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%%
% Regression of a, b is function of eta - 1

etas = linspace(2, 3, 11);
Weps = zeros(1, numel(epsilon));
aVec = zeros(1, numel(etas));
bVec = zeros(1, numel(etas));

for j=1:numel(etas)
    eta = etas(j);
    for i=1:numel(epsilon)
        er = epsilon(i);
        Weps(i) = potentialEnergyByRangeOfEta(er, eta, N, K, 1e-3);
    end
    
    [fitresult, gof] = createFit(epsilon, Weps);
    disp(gof);
    
    aVec(j) = fitresult.a;
    bVec(j) = fitresult.b;
end

figure; hold on;
plot(etas - 1, aVec, '--o', 'LineWidth', 2);
xlabel('$\eta - 1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$a(\eta-1)$', 'FontSize', 14, 'Interpreter', 'latex');

figure; hold on;
plot(etas - 1, bVec, '--o', 'LineWidth', 2);
xlabel('$\eta - 1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$b(\eta-1)$', 'FontSize', 14, 'Interpreter', 'latex');

function [fitresult, gof] = createFit(epsilon, Weps)
    %CREATEFIT(EPSILON,WEPS)
    %  Create a fit.
    %
    %  Data for 'untitled fit 1' fit:
    %      X Input : epsilon
    %      Y Output: Weps
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.
    %
    %  See also FIT, CFIT, SFIT.

    %  Auto-generated by MATLAB on 22-Jun-2021 12:28:48


    %% Fit: 'untitled fit 1'.
    [xData, yData] = prepareCurveData( epsilon, Weps );

    % Set up fittype and options.
    ft = fittype( 'a*(x-1)+b*(x-1)/(x+1)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.88716858638465 0.15193522204658];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end
