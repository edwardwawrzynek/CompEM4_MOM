function sigma = bistatic_echo_width_te_pec(phi, k, a)
    % truncation order of series
    N = 50; % 50 harmonics
    n = -N:1:N;

    % initialize
    sigma = zeros(size(phi));
    
    for p = 1:numel(phi)
        sigma(p) = 4/k * abs(sum(dbesselj(n, k*a) ./ dbesselh(n,2,k*a) .* exp(1j*n*phi(p)))).^2;
    end
end