% First-order derivative of the Modified Bessel Function Kn(x)

function deriv = besselkp(n, x)
    if (n == 0)
        deriv = -besselk(1, x);
    else
        deriv = -0.5 * (besselk(n-1, x) + besselk(n+1, x));
    end
end