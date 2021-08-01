% Potential energy on the point charge in the limit of a conducting
% cylinder.

function W = potentialEnergyOnPointChargeConducting(eta, N, K0, K, RelTol)
    C = 2 / pi;
    func = @(n, k, epsilonR, eta) potentialEnergyOnPointChargeConductingCoeff(n, k, epsilonR, eta);
    W = sumOfIntegralsSingle(C, func, 1, eta, N, K0, K, RelTol);
end