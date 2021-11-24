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
epsilon = 2;
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

epsilonVec = [4, 6, 8, 12];
resonances = {
    [3.376, 4, 4.638, 5.239, 5.84, 6.411, 6.982];
    [1.939, 2.405, 2.911, 3.411, 3.9, 4.38, 4.857, 5.324, 5.787];
    [1.678, 2.13, 2.568, 3, 3.421, 3.833, 4.243, 4.646, 5.044, 5.542];
    [1.412, 1.777, 2.134, 2.481, 2.824, 3.157, 3.488, 3.84, 4.212, 4.569, 4.921, 5.269, 5.557]
};
p1 = zeros(1, numel(epsilonVec));

disp('Resonance gof:')

figure; hold on;
for i=1:numel(epsilonVec)
    epsilon = epsilonVec(i);
    resonancesCurr = resonances{i};
    [fitresult, gof] = fit((1:numel(resonancesCurr)).', resonancesCurr.', 'poly1');
    p1(i) = fitresult.p1;
    disp(gof.rsquare);
    plot(resonancesCurr, '--o', 'LineWidth', 1, 'DisplayName', sprintf('$\\varepsilon=%d$', epsilon));
end

xlabel('$\#$Peak', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\omega R / c$', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontSize', 14, 'Interpreter', 'latex');

figure;
plot(epsilonVec, p1, '--o', 'LineWidth', 1);
xlabel('$\varepsilon$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$p_1$', 'FontSize', 14, 'Interpreter', 'latex');

[fitresult, gof] = fit(log10(epsilonVec.'), log10(p1.'), 'poly1');
slope = fitresult.p1;
disp('p1 gof:')
disp(gof.rsquare);