% Plot spectrum components

n = 1;
beta = 0.49;
omegas = linspace(0.01, 3, 200);
W = zeros(1, numel(omegas));
for i=1:numel(omegas)
    disp(omega);
            omega = omegas(i);
            kRange = linspace(1e-4, omega, M+1); kRange(end) = [];
            [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kRange, omega, x0, y0, z0, beta, epsilon);
            W(i) = omega .* trapz(kRange, 1./(omega^2 - kRange.^2) .* (abs(Bnk).^2 + abs(eta0Dnk).^2));
end

figure;
semilogy(omegas, W, 'LineWidth', 1); hold on;