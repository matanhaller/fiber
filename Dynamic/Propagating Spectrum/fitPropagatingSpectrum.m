% Regression on various propagating spectrum measurements.

% Clearing workspace
close all;
clear;
clc;

%% Finding exponential decay rate of spectrum w.r.t eta
epsilon = 4;
beta = 0.8;
gamma = (1 - beta^2)^(-0.5);
gb = gamma * beta;
M = 400;
omega = linspace(0.0101, 12.0001, M);
y0Vec = logspace(log10(1.05), 1, 40);
spec = zeros(numel(y0Vec), numel(omega));

fig = openfig(sprintf('spec_multiple_eta_eps_%d', epsilon), 'invisible');
dataObjs = findobj(fig, '-property', 'YData');

for i=1:numel(y0Vec)
    spec(i,:) = dataObjs(numel(y0Vec)+1-i).YData;
end

specl = log(spec);

omegaMin = 2.5;
omegaMax = 8;

y = specl(:, omega > omegaMin & omega < omegaMax);
x = omega(omega > omegaMin & omega < omegaMax);

p1 = zeros(1, numel(y0Vec));
p2 = zeros(1, numel(y0Vec));
for i=1:numel(p1)
    [fitresult, gof] = fit(x.', y(i,:).', 'poly1');
    p1(i) = fitresult.p1;
    p2(i) = fitresult.p2;
    disp(gof.rsquare);
end

figure;
plot(y0Vec - 1, p1, 'LineWidth', 1);
xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$p_1$', 'FontSize', 14, 'Interpreter', 'latex');

figure; hold on;
for i=1:numel(y0Vec)
    plot(omega, specl(i,:) - (p1(i)*omega+p2(i)), 'LineWidth', 1);
end
xlim([2, 8.5]);
ylim([-3, 5]);
xlabel('$\omega$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Subtracted Spectrum', 'FontSize', 14, 'Interpreter', 'latex');

x1 = y0Vec - 1;
y1 = p1;

[fitresult, gof] = fit(x1.', y1.', 'poly1');
disp(['p1 = ', num2str(fitresult.p1)]);
disp(['p2 = ', num2str(fitresult.p2)]);
disp(['rsquare = ', num2str(gof.rsquare)]);
disp(['p1*gb = ', num2str(fitresult.p1*gb)]);

%% Finding exponential decay rate of spectrum w.r.t beta
epsilon = 12;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = (0.1:0.1:4) * gbC; 
M = 400;
omega = linspace(0.0101, 12.0001, M);
betaVec = sqrt(gb.^2 ./ (gb.^2 + 1));
spec = zeros(numel(betaVec), numel(omega));

fig = openfig(sprintf('spec_multiple_beta_eps_%d', epsilon), 'invisible');
dataObjs = findobj(fig, '-property', 'YData');

for i=1:numel(betaVec)
    spec(i,:) = dataObjs(numel(betaVec)+1-i).YData;
end

specl = log(spec);

omegaMin = 4;
omegaMax = 8;

y = specl(:, omega > omegaMin & omega < omegaMax);
x = omega(omega > omegaMin & omega < omegaMax);

p1 = zeros(1, numel(betaVec));
for i=1:numel(p1)
    [fitresult, gof] = fit(x.', y(i,:).', 'poly1');
    p1(i) = fitresult.p1;
    disp(gof.rsquare);
end

xVec = 1./gb; xVec = xVec(8:end);

figure;
plot(xVec, p1(8:end), 'LineWidth', 1);
xlabel('$1/\gamma\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$p_1$', 'FontSize', 14, 'Interpreter', 'latex');

x1 = xVec;
y1 = p1(8:end);

[fitresult, gof] = fit(x1.', y1.', 'poly1');
disp(['p1 = ', num2str(fitresult.p1)]);
disp(['p2 = ', num2str(fitresult.p2)]);
disp(['rsquare = ', num2str(gof.rsquare)]);

%% Finding dependence of radiated energy w.r.t. gb
epsilon = 12;
betaC = 1 / sqrt(epsilon);
gbC = betaC / sqrt(1 - betaC^2);
gb = (0.1:0.1:4) * gbC; 
M = 400;
omega = linspace(0.0101, 12.0001, M);
betaVec = sqrt(gb.^2 ./ (gb.^2 + 1));
spec = zeros(numel(betaVec), numel(omega));

fig = openfig(sprintf('W_multiple_beta_eps_%d', epsilon), 'invisible');
dataObjs = findobj(fig, '-property', 'YData');

% Low gb
gbMax = 0.4*gbC;

Wl = log(dataObjs(1).YData); Wl = Wl(gb < gbMax);
xl = log(gb(gb < gbMax));

[fitresult, gof] = fit(xl.', Wl.', 'poly1');
disp(['p1 = ', num2str(fitresult.p1)]);
disp(['p2 = ', num2str(fitresult.p2)]);
disp(['rsquare = ', num2str(gof.rsquare)]);

% High gb
gbMin = 3*gbC;

W = dataObjs(1).YData; W = W(gb > gbMin);
x = gb(gb > gbMin);

[fitresult, gof] = fit(x.', W.', 'poly1');
disp(['p1 = ', num2str(fitresult.p1)]);
disp(['p2 = ', num2str(fitresult.p2)]);
disp(['rsquare = ', num2str(gof.rsquare)]);

%% Plotting radiated energy w.r.t gb
figure(6); hold on;

for epsilon=[2, 4, 12]
    betaC = 1 / sqrt(epsilon);
    gbC = betaC / sqrt(1 - betaC^2);
    gb = (0.1:0.1:4) * gbC;
    fig = openfig(sprintf('W_multiple_beta_eps_%d', epsilon), 'invisible');
    dataObjs = findobj(fig, '-property', 'YData');
    W = dataObjs(1).YData;
    
    figure(6);
    plot(gb, W, 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%d$', epsilon));
end

xlabel('$\gamma\beta$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\bar{W}$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

%% Finding resonance locations
close all;
clear;
clc;

% epsilonVec = [4, 6, 8, 12];
% resonances = {
%     [3.376, 4, 4.638, 5.239, 5.84, 6.411, 6.982];
%     [1.939, 2.405, 2.911, 3.411, 3.9, 4.38, 4.857, 5.324, 5.787];
%     [1.678, 2.13, 2.568, 3, 3.421, 3.833, 4.243, 4.646, 5.044, 5.542];
%     [1.412, 1.777, 2.134, 2.481, 2.824, 3.157, 3.488, 3.84, 4.212, 4.569, 4.921, 5.269, 5.557]
% };

M = 400;
omega = linspace(0.0101, 12.0001, M);

epsilonVec = linspace(2, 12, 40);
p1 = zeros(1, numel(epsilonVec));

fig = openfig('multiple_epsilon_eta_1.4_beta_0.8.fig', 'invisible');
dataObjs = findobj(fig, '-property', 'YData');

for i=1:numel(epsilonVec)
    spec(i,:) = dataObjs(numel(epsilonVec)+1-i).YData;
end

resonances = {
    [];
    [6.59, 7.31, 8.124, 8.844];
    [6.291, 7.042, 7.76, 8.484];
    [5.359, 6.08, 6.77, 7.463, 8.184];
    [5.149, 5.9, 6.56, 7.222, 7.883];
    [4.367, 5.028, 5.689, 6.32, 6.982, 7.613];
    [3.616, 4.247, 4.878, 5.509, 6.14, 6.77, 7.372];
    [3.556, 4.157, 4.758, 5.359, 5.96, 6.56, 7.16];
    [3.435, 4.037, 4.607, 5.209, 5.81, 6.38, 6.952];
%     [3.31, 3.947, 4.517, 5.088, 5.69, 6.2, 6.77, 7.31, 7.91, 8.514];
    [3.29, 3.857, 4.397, 4.968, 5.509, 6.05, 6.59, 7.13, 7.732, 8.304];
    [2.685, 3.165, 3.766, 4.307, 4.788, 5.359, 5.9, 6.411, 6.982, 7.553, 8.094];
    [2.624, 3.135, 3.616, 4.217, 4.728, 5.239, 5.78, 6.261, 6.831, 7.372, 7.913];
    [2.054, 2.564, 3.075, 3.586, 4.127, 4.638, 5.149, 5.629, 6.14, 6.681, 7.222, 7.732];
    [2, 2.504, 3.015, 3.496, 4.037, 4.547, 5.028, 5.509, 6, 6.53, 7.072, 7.583];
    [1.963, 2.444, 2.955, 3.465, 3.946, 4.457, 4.938, 5.389, 5.87, 6.41, 6.922, 7.43];
    [1.482, 1.933, 2.414, 2.895, 3.406, 3.887, 4.367, 4.818, 5.359, 5.81, 6.291, 6.771, 7.282, 7.67, 8.154, 8.664];
    [1.452, 1.873, 2.384, 2.83, 3.316, 3.796, 4.277, 7.728, 5.209, 5.629, 6.17, 6.65, 7.132, 7.522, 8, 8.484];
    [1.422, 1.843, 2.323, 2.804, 3.285, 3.736, 4.187, 4.668, 5.089, 5.569, 6.05, 6.53, 6.95, 7.34, 7.85, 8.334];
    [0.972, 1.39, 1.813, 2.294, 2.745, 3.23, 3.68, 4.127, 4.608, 5, 5.45, 5.96, 6.41, 6.892, 7.19, 7.703, 8.184];
    [0.972, 1.36, 1.78, 2.264, 2.715, 3.16, 3.616, 4.067, 4.49, 4.9, 5.389, 5.84, 6.32, 6.77, 7.102, 7.583, 8.064];
    [0.911, 1.33, 1.753, 2.203, 2.654, 3.105, 3.556, 4, 4.36, 5.3, 5.75, 6.23, 6.44, 6.68, 6.982, 7.462, 7.91];
    [0.911, 1.33, 1.75, 2.17, 2.62, 3.075, 3.496, 3.912, 4.337, 4.758, 5.209, 5.66, 6.11, 6.411, 6.56, 6.89, 7.34, 7.793];
    [0.912, 1.302, 1.723, 2.143, 2.594, 3.015, 3.466, 3.85, 4.217, 4.668, 5.119, 5.569, 6.02, 6.32, 6.47, 6.77, 7.222, 7.67];
    [0.88, 1.27, 1.69, 2.11, 2.56, 2.985, 3.406, 3.82, 4.157, 4.608, 5.058, 5.509, 5.93, 6.35, 6.68, 7.132, 7.583, 7.85, 8.274, 8.725];
    [0.88, 1.27, 1.66, 2.083, 2.504, 2.955, 3.346, 3.76, 4.157, 4.547, 5, 5.419, 5.839, 6.11, 6.56, 7.01, 7.462, 8.154, 8.604];
    [0.851, 1.24, 1.633, 2.05, 2.47, 2.895, 3.19, 3.61, 4.037, 4.488 4.908, 5.359, 5.78, 6.05, 6.471, 6.9212, 7.372, 8.03, 8.484];
    [0.851, 1.212, 1.632, 2.023, 2.444, 2.864, 3.165, 3.556, 4.037, 4.427, 4.848, 5.268, 5.479, 5.96, 6.381, 6.83, 7.252, 7.91, 8.36];
    [0.821, 1.212, 1.603, 2.023, 2.313, 2.83, 3.105, 3.526, 3.945, 4.367, 4.79, 5.209, 5.449, 5.867, 6.32, 6.74, 7.16, 7.82, 8.244];
    [0.82, 1.18, 1.572, 2, 2.38, 2.83, 3.07, 3.465, 3.886, 4.307, 4.727, 5.149, 5.33, 5.81, 6.23, 6.65, 7.072, 7.703, 8.15];
    [0.79, 1.18, 1.57, 1.96, 2.35, 2.74, 3.04, 3.436, 3.857, 4.247, 4.668, 5.059, 5.329, 5.72, 6.14, 6.56, 6.982, 7.22, 7.61, 8.03];
    [0.79, 1.15, 1.54, 1.93, 2.32, 2.72, 2.98, 3.37, 3.8, 4.22, 4.6, 5, 5.23, 5.66, 6.08, 6.5, 6.89, 7.52, 7.94];
    [0.791, 1.122, 1.512, 1.903, 2.294, 2.564, 2.955, 3.34, 3.73, 4.157, 4.548, 4.94, 5.179, 5.6, 6.02, 6.412, 6.8, 7.042, 7.43, 7.85];
    [0.76, 1.12, 1.51, 1.9, 2.294, 2.53, 2.92, 3.31, 3.73, 4.097, 4.51, 4.878, 5.12, 5.54, 5.93, 6.35, 6.86, 7.34, 7.76, 8.36, 8.75];
    [0.76, 1.12, 1.48, 1.87, 2.054, 2.26, 2.89, 3.255, 3.67, 4.067, 4.457, 4.667, 5.058, 5.79, 5.869, 6.26, 6.65, 6.86, 7.28, 7.67, 8.27, 8.664];
%     [0.761, 1.12, 1.48, 1.84, 2.23, 2.864, 3.225, 3.616, 4.01, 4.397, 4.608, 5, 5.419, 5.81, 6.2, 6.77, 7.19, 7.58, 8.18, 8.57];
    [0.731, 1.092, 1.452, 1.84, 2.204, 2.44, 2.83, 3.19, 3.58, 3.97, 4.367, 4.578, 4.938, 5.36, 5.75, 6.14, 6.65, 7.1, 7.522, 8.06, 8.484];
%     [0.731, 1.09, 1.45, 1.813, 2.174, 2.414, 2.774, 3.556, 3.947, 4.517, 4.908, 5.3, 5.69, 6.08, 6.65, 7.04, 7.43, 7.97, 8.394];
    [0.731, 1.062, 1.422, 1.813, 2.173, 2.384, 2.774, 3.135, 3.525, 3.887, 4.277, 4.45, 4.848, 5.239, 5.629, 6, 6.56, 6.951, 7.34, 7.91, 8.3];
    [0.701, 1.06, 1.4, 1.783, 2.14, 2.354, 2.744, 3.105, 3.466, 3.857, 4.217, 4.427, 4.788, 5.179, 5.569, 5.93, 6.5, 6.89, 7.28, 7.82, 8.214];
};

for i=1:numel(epsilonVec)-3
    resonances{i} = sort(resonances{i});
    resonances{i}(resonances{i} > 5) = [];
end

for i=1:numel(epsilonVec)-3
    epsilon = epsilonVec(i);
    resonancesCurr = resonances{i};
    if numel(resonancesCurr) > 2
        [fitresult, gof] = fit((1:numel(resonancesCurr)).', resonancesCurr.', 'poly1');
        p1(i) = fitresult.p1;
        disp(gof.rsquare);
    end
end

figure; hold on;
plot(epsilonVec, p1, 'o', 'LineWidth', 1);
xlabel('$\varepsilon$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$p_1$', 'FontSize', 14, 'Interpreter', 'latex');
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');

epsilonVecFit = epsilonVec(p1 > 0);
p1Fit = p1(p1 > 0);

[fitresult, gof] = fit(log10((epsilonVecFit).'), log10(p1Fit.'), 'poly1');
slope = fitresult.p1;
disp('p1 gof:')
disp(gof.rsquare);

x = linspace(2, 12, 1e3);
plot(x, 10^fitresult.p2 * x.^slope, 'LineWidth', 1);