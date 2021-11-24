% Finding cutoff frequency of each eigenmode of the dielectric cylinder
% according to the conditions in 'The Essence of Dielectric Waveguides' (p.
% 146)

function wco = findCutoffFreq(n, s, epsilon)
    m = ceil(s / 2); % Eigenmode index
    t = (mod(s, 2) == 0); % Eigenmode type (TE/TM or HE/EH)
    
    if n == 0
        wco = fzero(@(z) besselj(n,z), besseljZeroGuess(n,m)) / sqrt(epsilon - 1);
    
    elseif n == 1
        if t == 0
            wco = fzero(@(z) besselj(n,z), besseljZeroGuess(n,m-1)) / sqrt(epsilon - 1);
        else
            wco = fzero(@(z) besselj(n,z), besseljZeroGuess(n,m)) / sqrt(epsilon - 1);
        end
    
    else
        if t == 0
            wco = fzero(@(z) z/(n-1).*besselj(n,z)-(epsilon+1).*besselj(n-1,z),...
                        besseljZeroGuess(n-1,m)) / sqrt(epsilon - 1);
        else
            wco = fzero(@(z) besselj(n,z), besseljZeroGuess(n,m)) / sqrt(epsilon - 1);
        end
    end
    
    wco = abs(wco); % Fixing negative values
end
    
% Guess of Jn(x) m-th zero
% (Credit to https://www.mathworks.com/matlabcentral/answers/230834-zeros-of-bessel-functions)
function guess = besseljZeroGuess(n, m)
    guess = 2.5505 + 1.2474 * n + (m - 1) * pi;
end