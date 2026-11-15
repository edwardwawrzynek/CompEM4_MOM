function sigma = bistatic_echo_width_tm_eps(phi, k, a, ep_r)
    % truncation order of series
    N = 50; % 50 harmonics
    n = -N:1:N;

    % initialize
    sigma = zeros(size(phi));
    
    kd = sqrt(ep_r)*k;
    
    for p = 1:numel(phi)
        alpha_n = (sqrt(ep_r).*besselj(n, k*a).*dbesselj(n, kd*a) - dbesselj(n, k*a).*besselj(n, kd*a)) ./ (sqrt(ep_r).*besselh(n,2,k*a).*dbesselj(n,kd*a) - dbesselh(n,2,k*a).*besselj(n,kd*a));   % CHANGED: removed (-1j).^(-n)
        sigma(p) = 4/k * abs(sum(alpha_n .* exp(1j*n*phi(p)))).^2;
    end
end