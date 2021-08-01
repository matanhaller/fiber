% First-order derivative of the Bessel Function Jn(x)

function deriv = besseljp(n, x)
    if (n == 0)
        deriv = -besselj(1, x);
    else
        deriv = 0.5 * (besselj(n-1, x) - besselj(n+1, x));
    end
end