% Template for expressions consisting of an infinite sum of integrals (of
% 3 arguments)

function res = sumOfIntegralsTriple(C, coeffFunc, epsilonR, eta, x, y, z, N, K, RelTol)
    res = 0;
    for n=0:N
        func = @(k) coeffFunc(n, k, epsilonR, eta, x, y, z);
        intOnK = integral(func, 1e-6, K, 'RelTol', RelTol, 'ArrayValued', true);
        if (n == 0)
            res = intOnK;
        else
            res = res + 2 * intOnK;
        end
    end
    res = res .* C;
end