% Calculates force according to eta's order of magnitude: if eta~1,
% calculate for Large N and large K0. If eta>>1, calculate for small N and
% small K0;

function phi = forceByRangeOfEta(epsilonR, eta, N, K, RelTol)
    etaLo = eta(eta < 1.1);
    etaHi = eta(eta >= 1.1);
    
    Nhi = 20;
    
    phi = [forceOnPointCharge(epsilonR, etaLo, N, 1e-2, K, RelTol), ...
           forceOnPointCharge(epsilonR, etaHi, min(Nhi, N), 0, K, RelTol)];
end