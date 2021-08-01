% Secondary potential induced by the fiber

function phi = secondaryPotential(epsilonR, eta, x, y, z, N, K, RelTol)
    C = -(2 / pi);
    func = @(n, k, epsilonR, eta, x, y, z) secondaryPotentialCoeff(n, k, epsilonR, eta, x, y, z);
    phi = sumOfIntegralsTriple(C, func, epsilonR, eta, x, y, z, N, K, RelTol);
end