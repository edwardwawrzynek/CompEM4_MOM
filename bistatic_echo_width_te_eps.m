function sigma = bistatic_echo_width_te_eps(phi, k, a, ep_r)
    % truncation order of series
    N = 50; % 50 harmonics
    n = -N:1:N;

    % initialize
    sigma = zeros(size(phi));
    
    kd = sqrt(ep_r)*k;
    
    alpha_n = (sqrt(ep_r).*dbesselj(n, k*a).*besselj(n, kd*a) - besselj(n, k*a).*dbesselj(n, kd*a)) ./ (sqrt(ep_r).*dbesselh(n,2,k*a).*besselj(n,kd*a) - besselh(n,2,k*a).*dbesselj(n,kd*a));   % CHANGED: removed (-1j).^(-n); moved out of loop (independent of phi)
    
    for p = 1:numel(phi)
        sigma(p) = 4/k * abs(sum(alpha_n .* exp(1j*n*phi(p)))).^2;
    end
end