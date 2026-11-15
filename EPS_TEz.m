clc; clear; close all;

eps0 = 8.85 * 1e-12;
mu0 = 4 * pi * 1e-7;
c0 = 1 / sqrt(eps0 * mu0);
f = 10 * 1e9;
lambda = c0/f;
k = 2*pi/lambda;
geo_a = 1.5 * 1e-2;
eps_r = 9;

% initialize the problem
a = CFIE_MOMProblem(k);

% create the geometry
a = a.makeCircle(geo_a, 100);

% boundary condition is not needed for the solution domain (yay!)

% assign incidence
a = a.constructIncidence(0, 1);

% assemble LHS
a = a.tez_eps_sys_matrix(eps_r);

% solve the method-of-moments system
a = a.solveI();

% construct the field
step_num = 100;
Xpoints = linspace(-5*geo_a, 5*geo_a, step_num);
Ypoints = linspace(-5*geo_a, 5*geo_a, step_num);
[Hinc, Hscat, Htot] = a.evaluateEPSFields(Xpoints, Ypoints, 2);

% % plot the mesh
% figure();
% a.plotMesh();
% set(gcf, 'Color', 'w');
% set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
% box on; grid on;
% axis equal tight;
% title("Mesh")
% xlabel("x [m]")
% ylabel("y [m]")

% plot the field
fig = figure('Color','White', 'Units','normalized', 'Position', [0.02 0.04 0.96 0.5]);
t = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');

nexttile; hold on;
imagesc(Xpoints, Ypoints, real(Hinc));
theta = linspace(0, 2*pi, 100);
plot(geo_a*cos(theta), geo_a*sin(theta), 'k', 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
set(gca, 'YDir', 'normal', 'FontName', 'times', 'FontSize', 18);
colormap inferno;
colorbar;
clim([-2,2]);
axis equal tight;
title("Incident Field [V/m]");
xlabel("x [m]")
ylabel("y [m]")

nexttile; hold on;
imagesc(Xpoints, Ypoints, real(Hscat));
theta = linspace(0, 2*pi, 100);
plot(geo_a*cos(theta), geo_a*sin(theta), 'k', 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
set(gca, 'YDir', 'normal', 'FontName', 'times', 'FontSize', 18);
colormap inferno;
colorbar;
clim([-2,2]);
axis equal tight;
title("Scattered Field [V/m]");
xlabel("x [m]")
ylabel("y [m]")

nexttile; hold on;
imagesc(Xpoints, Ypoints, real(Htot));
theta = linspace(0, 2*pi, 100);
plot(geo_a*cos(theta), geo_a*sin(theta), 'k', 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
set(gca, 'YDir', 'normal', 'FontName', 'times', 'FontSize', 18);
colormap inferno;
colorbar;
clim([-2,2]);
axis equal tight;
title("Total Field [V/m]");
xlabel("x [m]")
ylabel("y [m]")

%%
%%%%%%%%%%%%%%%%%%
%%% Echo Width %%%
%%%%%%%%%%%%%%%%%%

% sweep the angle
angle = linspace(0, 2*pi, 360);

[X, Y] = meshgrid(Xpoints, Ypoints);
X = X(:); Y = Y(:);
sigma_simulated  = simulated_echo_width(X, Y, k, Htot, 0.05, 360, 50);
sigma_analytical = bistatic_echo_width_te_eps(angle, k, geo_a, eps_r);

% plot
figure();
plot(angle*180/pi, sigma_simulated, 'LineWidth', 2); hold on;
plot(angle*180/pi, sigma_analytical, ":", 'LineWidth', 2);
set(gca,'FontName', 'Times', 'FontSize', 18);
set(gcf, 'Color', 'w')
legend("Simulated", "Analytic");
title("Bistatic Echowidth; TE_z Dielectric")
xlabel("\phi [deg]");
ylabel("\sigma(\phi)");
grid on; box on;