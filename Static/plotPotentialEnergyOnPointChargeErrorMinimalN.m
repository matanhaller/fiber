% Finding the minimal order which yields a relative error lower than a
% given threshold.

close all;
clear;
clc;

epsilonR = [1.2, 2, 4, 12];
Ninf = 80;
K = 500;

etas = logspace(log10(1.02), 2, 100);
thresholds = [5e-2, 1e-2, 0.1e-2];

figure;
tiledlayout(2, 2);

for er=epsilonR
    minN = Ninf * ones(numel(thresholds), numel(etas));
    for i=1:numel(etas)
        eta = etas(i);
        phiInf = potentialEnergyByRangeOfEta(er, eta, Ninf, K, 1e-3);
        n = 0;
        for j=1:numel(thresholds)
            t = thresholds(j);
            err = 1;
            while err >= t
                disp(n);
                if n > Ninf
                    break
                end
                phiTrunc = potentialEnergyByRangeOfEta(er, eta, n, K, 1e-3);
                err = relError(phiInf, phiTrunc);
                n = n + 1;
            end
            n = n - 1;
            minN(j, i) = n;
        end
    end
    
    nexttile;
    for i=1:numel(thresholds)
        loglog(etas - 1, minN(i,:), 'LineWidth', 2, 'DisplayName', sprintf('$\\mathrm{Error} < %.1f$\\%%', thresholds(i) * 100)); hold on;
    end
    title(sprintf('$\\varepsilon=%.1f$', er), 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('$\eta-1$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Minimal Order', 'FontSize', 14, 'Interpreter', 'latex');
    legend('FontSize', 14, 'Interpreter', 'latex');
end