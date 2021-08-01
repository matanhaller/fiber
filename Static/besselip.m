% First-order derivative of the Modified Bessel Function In(x)

function deriv = besselip(n, x)
    if (n == 0)
        deriv = besseli(1, x);
    else
        deriv = 0.5 * (besseli(n-1, x) + besseli(n+1, x));
    end
end