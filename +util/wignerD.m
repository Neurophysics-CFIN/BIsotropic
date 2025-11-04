function D = wignerD(l, m1, m2, alpha, beta, gamma)
    % Validate inputs
    if abs(m1) > l || abs(m2) > l
        error('Invalid m1 or m2 for given l');
    end
    
    % Calculate the Wigner small-d matrix element
    d = wigner_small_d(l, m1, m2, beta);
    
    % Calculate the Wigner D-matrix element
    D = exp(-1i * m1 * alpha) .* d .* exp(-1i * m2 * gamma);
end

function d = wigner_small_d(l, m1, m2, beta)
    % This function calculates the Wigner small-d matrix element d^l_{m1,m2}(beta)
    d = 0;
    % Calculate the prefactor
    prefactor = sqrt(factorial(l + m1) * factorial(l - m1) * factorial(l + m2) * factorial(l - m2));
    for k = max(0, m2 - m1) : min(l + m2, l - m1)
            % Calculate the individual terms in the sum
            numerator = (-1)^(m1 - m2 + k).*(cos(beta / 2).^(2 * l + m2 - m1 - 2 * k)) .* (sin(beta / 2).^(m1 - m2 + 2 * k));
            denominator = factorial(k) * factorial(m1 - m2 + k) * factorial(l + m2 - k) * factorial(l - m1 - k);
            d = d + (numerator / denominator);
        % end
    end
    d = prefactor * d;
end
