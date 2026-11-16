clc; clear; close all;

sweep_var = 3;

if sweep_var == 1
    % number of segments in the cylinder mesh
    sweep = [50, 100, 200, 400, 800];
end

if sweep_var == 2
    % solution domain size relative to the wavelength
    sweep = [5, 10, 15, 20, 25];
end

if sweep_var == 3
    % number of steps to which evaluate the field in x, y directions
    sweep = [50, 100, 200, 400, 800];
end

% initialize
RMSE_PEC_TMz = zeros(1, numel(sweep));
RMSE_PEC_TEz = zeros(1, numel(sweep));
RMSE_EPS_TMz = zeros(1, numel(sweep));
RMSE_EPS_TEz = zeros(1, numel(sweep));

for idx = 1:numel(sweep)

    if sweep_var == 1
        number_of_segments = sweep(idx);
        domain_bound = 10;
        step_num = 200;
    end

    if sweep_var == 2
        number_of_segments = 200;
        domain_bound = sweep(idx);
        step_num = 200;
    end

    if sweep_var == 3
        number_of_segments = 200;
        domain_bound = 10;
        step_num = sweep(idx);
    end

    eps0 = 8.85 * 1e-12;
    mu0 = 4 * pi * 1e-7;
    c0 = 1 / sqrt(eps0 * mu0);
    f = 10 * 1e9;
    lambda = c0/f;
    k = 2*pi/lambda;
    geo_a = 1.5 * 1e-2;
    eps_r = 9;

    %%% Run PEC TMz
    % initialize the problem
    a = MOMProblem(k);

    % create the geometry
    a = a.makeCircle(geo_a, number_of_segments);

    % boundary condition is not needed for the solution domain (yay!)

    % assign incidence
    a = a.constructIncidence(0, 1);

    % assemble LHS
    a = a.tmz_pec_sys_matrix();

    % solve the method-of-moments system
    a = a.solveI();

    % construct the field
    Xpoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    Ypoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    [Einc_PEC, Escat_PEC, Etot_PEC] = a.evaluatePECFields(Xpoints, Ypoints, 1);

    %%% Run PEC TEz
    % initialize the problem
    b = MOMProblem(k);

    % create the geometry
    b = b.makeCircle(geo_a, 100);

    % boundary condition is not needed for the solution domain (yay!)

    % assign incidence
    b = b.constructIncidence(0, 1);

    % assemble LHS
    b = b.tez_pec_sys_matrix();

    % solve the method-of-moments system
    b = b.solveI();

    % construct the field
    Xpoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    Ypoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    [Hinc_PEC, Hscat_PEC, Htot_PEC] = b.evaluatePECFields(Xpoints, Ypoints, 2);

    %%% Run Dielectric TMz
    % initialize the problem
    c = CFIE_MOMProblem(k);

    % create the geometry
    c = c.makeCircle(geo_a, 100);

    % boundary condition is not needed for the solution domain (yay!)

    % assign incidence
    c = c.constructIncidence(0, 1);

    % assemble LHS
    c = c.tmz_eps_sys_matrix(eps_r);

    % solve the method-of-moments system
    c = c.solveI();

    % construct the field
    Xpoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    Ypoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    [Einc_EPS, Escat_EPS, Etot_EPS] = c.evaluateEPSFields(Xpoints, Ypoints, 1);

    %%% Run Dielectric TEz
    % initialize the problem
    d = CFIE_MOMProblem(k);

    % create the geometry
    d = d.makeCircle(geo_a, number_of_segments);

    % boundary condition is not needed for the solution domain (yay!)

    % assign incidence
    d = d.constructIncidence(0, 1);

    % assemble LHS
    d = d.tez_eps_sys_matrix(eps_r);

    % solve the method-of-moments system
    d = d.solveI();

    % construct the field
    Xpoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    Ypoints = linspace(-domain_bound*lambda, domain_bound*lambda, step_num);
    [Hinc_EPS, Hscat_EPS, Htot_EPS] = d.evaluateEPSFields(Xpoints, Ypoints, 2);

    %%% Save the RMSE
    angle = linspace(0, 2*pi, 360);
    [X, Y] = meshgrid(Xpoints, Ypoints);
    X = X(:); Y = Y(:);

    % PEC TMz
    sigma_simulated  = simulated_echo_width(X, Y, k, Etot_PEC, 0.05, 360, 50);
    sigma_analytical = bistatic_echo_width_tm_pec(angle, k, geo_a);
    RMSE_PEC_TMz(idx) = ...
        sqrt(sum((sigma_simulated(:) - sigma_analytical(:)).^2)) ...
        ./ sqrt(sum(sigma_analytical(:).^2));

    % PEC TEz
    sigma_simulated  = simulated_echo_width(X, Y, k, Htot_PEC, 0.05, 360, 50);
    sigma_analytical = bistatic_echo_width_te_pec(angle, k, geo_a);
    RMSE_PEC_TEz(idx) = ...
        sqrt(sum((sigma_simulated(:) - sigma_analytical(:)).^2)) ...
        ./ sqrt(sum(sigma_analytical(:).^2));

    % Dielectric TMz
    sigma_simulated = simulated_echo_width(X, Y, k, Etot_EPS, 0.05, 360, 50);
    sigma_analytical = bistatic_echo_width_tm_eps(angle, k, geo_a, eps_r);
    RMSE_EPS_TMz(idx) = ...
        norm(sigma_simulated(:) - sigma_analytical(:)) / norm(sigma_analytical(:));

    % Dielectric TEz
    sigma_simulated = simulated_echo_width(X, Y, k, Htot_EPS, 0.05, 360, 50);
    sigma_analytical = bistatic_echo_width_te_eps(angle, k, geo_a, eps_r);
    RMSE_EPS_TEz(idx) = ...
        norm(sigma_simulated(:) - sigma_analytical(:)) / norm(sigma_analytical(:));

end

%% Plotting RMSE Sweeps

if sweep_var == 1
    figure(); hold on;
    mesh_size = sweep/lambda;
    plot(mesh_size, RMSE_PEC_TMz, 'LineWidth', 2);
    plot(mesh_size, RMSE_PEC_TEz, ":", 'LineWidth', 2);
    plot(mesh_size, RMSE_EPS_TMz, "--", 'LineWidth', 2);
    plot(mesh_size, RMSE_EPS_TEz, "-.", 'LineWidth', 2);
    set(gca,'FontName', 'Times', 'FontSize', 18);
    set(gca,'XScale', 'log', 'YScale', 'log');
    ylim([1e-2,1e1]);
    set(gcf, 'Color', 'w');
    legend(["TM_z PEC", "TE_z PEC", "TM_z Dielectric", "TE_z Dielectric"], 'Location', 'northwest');
    title("Bistatic Echowidth RMSE")
    xlabel("Mesh Size /\lambda");
    ylabel("RMSE");
    grid on; box on;
end

if sweep_var == 2
    figure(); hold on;
    plot(sweep, RMSE_PEC_TMz, 'LineWidth', 2);
    plot(sweep, RMSE_PEC_TEz, ":", 'LineWidth', 2);
    plot(sweep, RMSE_EPS_TMz, "--", 'LineWidth', 2);
    plot(sweep, RMSE_EPS_TEz, "-.", 'LineWidth', 2);
    set(gca,'FontName', 'Times', 'FontSize', 18);
    set(gca,'XScale', 'log', 'YScale', 'log');
    ylim([1e-2,1e1]);
    set(gcf, 'Color', 'w');
    legend(["TM_z PEC", "TE_z PEC", "TM_z Dielectric", "TE_z Dielectric"], 'Location', 'northwest');
    title("Bistatic Echowidth RMSE")
    xlabel("Solution Domain \lambda");
    ylabel("RMSE");
    grid on; box on;
end

if sweep_var == 3
    figure(); hold on;
    plot(sweep, RMSE_PEC_TMz, 'LineWidth', 2);
    plot(sweep, RMSE_PEC_TEz, ":", 'LineWidth', 2);
    plot(sweep, RMSE_EPS_TMz, "--", 'LineWidth', 2);
    plot(sweep, RMSE_EPS_TEz, "-.", 'LineWidth', 2);
    set(gca,'FontName', 'Times', 'FontSize', 18);
    set(gca,'XScale', 'log', 'YScale', 'log');
    ylim([1e-2,1e1]);
    set(gcf, 'Color', 'w');
    legend(["TM_z PEC", "TE_z PEC", "TM_z Dielectric", "TE_z Dielectric"], 'Location', 'northwest');
    title("Bistatic Echowidth RMSE")
    xlabel("Number of Grid Points");
    ylabel("RMSE");
    grid on; box on;
end