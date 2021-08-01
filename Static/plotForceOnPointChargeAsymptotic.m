% Asymptotic analysis of the force on the point charge due to
% the fiber

close all;
clear;
clc;

% Plotting for different values of epsilonR
epsilonR = 4;
eta = logspace(log10(1.02), log10(1e3), 100);
N = 80;
K = 500;

figure;

F = forceByRangeOfEta(epsilonR, eta, N, K, 1e-3);
loglog(eta - 1, -F, 'LineWidth', 2, 'DisplayName', 'Numeric');
hold on;

% for n=0:1
%     Ftrunc = forceByRangeOfEta(epsilonR, eta, n, K, 1e-3);
%     loglog(eta - 1, -Ftrunc, 'LineWidth', 2, 'DisplayName', sprintf('$N=%d$', n));
% end

% Dielectric half-plane force (eta-->1)
F0 = forceDielectricHalfPlane(epsilonR, eta);
loglog(eta - 1, F0, '--', 'LineWidth', 2, 'DisplayName', 'Asymptotic $\propto (\eta -1)^{-2}$');

% Induced dipole with ~1/(eta-1) dipole moment force (eta-->inf)
a = 1.5 * 0.589;
b = 1.5 * 0.098;
Finf = ((epsilonR - 1) / (epsilonR + 1) * a + (epsilonR - 1) * b) * (eta - 1).^(-4);
loglog(eta - 1, Finf, '--', 'LineWidth', 2, 'DisplayName', 'Asymptotic $\propto  (\eta -1)^{-4}$');

legend('FontSize', 14, 'Interpreter', 'latex');
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$F$', 'FontSize', 14, 'Interpreter', 'latex');

p = 0.5;
Ftrans = (F0.^(-p) + Finf.^(-p)).^(-1/p);
loglog(eta - 1, Ftrans, '--', 'LineWidth', 2, 'DisplayName', 'Analytic');

% Plotting relative error
epsilonR = [1.2, 2, 4, 12];
eta = logspace(log10(1.02), log10(1e3), 100);

figure;  
loglog(eta - 1, relError(-F, Ftrans), 'LineWidth', 2);
xlabel('$\eta - 1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Error', 'FontSize', 14, 'Interpreter', 'latex');

% Plotting for variable permittivity
epsilon = logspace(log10(1.01), log10(12), 40);
eta = 2.5;
Feps = zeros(1, numel(epsilon));

for i=1:numel(epsilon)
    er = epsilon(i);
    Feps(i) = forceByRangeOfEta(er, eta, N, K, 1e-3);
end

[fitresult, gof] = createFit(epsilon, -Feps);
disp(gof);

a = fitresult.a;
b = fitresult.b;

Ffit = a * (epsilon - 1) + b * (epsilon - 1) ./ (epsilon + 1);

figure; hold on;
plot(epsilon, -Feps, 'LineWidth', 2, 'DisplayName', 'Numeric');
plot(epsilon, Ffit, '--', 'LineWidth', 2, 'DisplayName', 'Analytic');
xlabel('$\varepsilon$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{F}$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%%
% Regression of a, b is function of eta - 1

etas = linspace(2, 3, 11);
Feps = zeros(1, numel(epsilon));
aVec = zeros(1, numel(etas));
bVec = zeros(1, numel(etas));

for j=1:numel(etas)
    eta = etas(j);
    for i=1:numel(epsilon)
        er = epsilon(i);
        Feps(i) = forceByRangeOfEta(er, eta, N, K, 1e-3);
    end
    
    [fitresult, gof] = createFit(epsilon, Feps);
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

function [fitresult, gof] = createFit(epsilon, Feps)
    %CREATEFIT(EPSILON,FEPS)
    %  Create a fit.
    %
    %  Data for 'untitled fit 1' fit:
    %      X Input : epsilon
    %      Y Output: Feps
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.
    %
    %  See also FIT, CFIT, SFIT.

    %  Auto-generated by MATLAB on 22-Jun-2021 12:28:48


    %% Fit: 'untitled fit 1'.
    [xData, yData] = prepareCurveData( epsilon, Feps );

    % Set up fittype and options.
    ft = fittype( 'a*(x-1)+b*(x-1)/(x+1)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.88716858638465 0.15193522204658];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end

