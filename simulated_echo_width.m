function sigma_simulated = simulated_echo_width(X, Y, k, Ftot, R_samp, Npts, Nmax)

    %%% Setup
    Einc = exp(-1j*k*X); % incident field
    Escat = Ftot(:) - Einc; % scattered field

    %%% Sample scattered field on circle 
    % points on the circle:
    angle = linspace(0, 2*pi, Npts+1);
    angle = angle(1:end-1); % exclude 2pi
    x = R_samp*cos(angle); 
    y = R_samp*sin(angle);

    % interpolate real and imag components
    Fr = scatteredInterpolant(X, Y, real(Escat), 'linear', 'none');
    Fi = scatteredInterpolant(X, Y, imag(Escat), 'linear', 'none');

    % apply over the circle and recombine
    Es_c = Fr(x, y) + 1j*Fi(x, y);

    % azimuthal harmonics
    a_n = fftshift(fft(Es_c)) / Npts;
    n = -Npts/2:Npts/2-1;

    % truncate to Nmax
    a_n = a_n(abs(n) <= Nmax);
    n = n(abs(n) <= Nmax);

    % propagate to far field
    c_n = a_n ./ besselh(n, 2, k*R_samp);

    % initialize
    sigma_simulated = zeros(size(angle));

    % calculate bistatic echo width
    for p = 1:numel(angle)
        sigma_simulated(p) = 4/k * abs(sum((1j).^n .* c_n .* exp(1j*n*angle(p))))^2;
    end
end