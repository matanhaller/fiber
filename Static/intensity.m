% Intensity of the E-field
function int = intensity(epsilonR, eta, x, y, z, N, K, RelTol)
    C = 2 / pi;
    func_1 = @(n, k, epsilonR, eta, x, y, z) intensityCoeffRadial(n, k, epsilonR, eta, x, y, z);
    radial = (sumOfIntegralsTriple(C, func_1, epsilonR, eta, x, y, z, N, K, RelTol)).^2;
    func_2 = @(n, k, epsilonR, eta, x, y, z) intensityCoeffAzimuthal(n, k, epsilonR, eta, x, y, z);
    azimuthal = (sumOfIntegralsTriple(C, func_2, epsilonR, eta, x, y, z, N, K, RelTol)).^2;
    func_3 = @(n, k, epsilonR, eta, x, y, z) intensityCoeffZ(n, k, epsilonR, eta, x, y, z);
    longitudenal = (sumOfIntegralsTriple(C, func_3, epsilonR, eta, x, y, z, N, K, RelTol)).^2;
    int = radial + azimuthal + longitudenal;
end
    