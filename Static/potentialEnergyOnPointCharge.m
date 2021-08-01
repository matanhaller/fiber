% Potential energy on the point charge

function W = potentialEnergyOnPointCharge(epsilonR, eta, N, K0, K, RelTol)
    C = 2 / pi;
    func = @(n, k, epsilonR, eta) potentialEnergyOnPointChargeCoeff(n, k, epsilonR, eta);
    W = sumOfIntegralsSingle(C, func, epsilonR, eta, N, K0, K, RelTol);
end