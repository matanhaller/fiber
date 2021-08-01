% Template for expressions consisting of an infinite sum of integrals (of
% single argument)

function res = sumOfIntegralsSingle(C, coeffFunc, epsilonR, eta, N, K0, K, RelTol)
    res = 0;
    for n=0:N
        func = @(k) coeffFunc(n, k, epsilonR, eta);
        intOnK = integral(func, K0, K, 'RelTol', RelTol, 'ArrayValued', true);
        if (n == 0)
            res = intOnK;
        else
            res = res + 2 * intOnK;
        end
    end
    res = res .* C;
end