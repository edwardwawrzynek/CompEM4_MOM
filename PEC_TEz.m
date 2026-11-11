clc; clear; close all;

eps0 = 8.85 * 1e-12;
mu0 = 4 * pi * 1e-7;
c0 = 1 / sqrt(eps0 * mu0);
f = 9 * 1e9;
lambda = c0/f;
k = 2*pi/lambda;
geo_a = 0.5 * lambda;

% initialize the problem
a = MOMProblem(k);

% create the geometry
a = a.makeCircle(geo_a, 100);

% boundary condition is not needed for the solution domain (yay!)

% assign incidence
a = a.constructIncidence(0, 1);

% assemble LHS
a = a.tez_pec_sys_matrix();

% solve the method-of-moments system
a = a.solveI();

% construct the field
step_size = 100;
Xpoints = linspace(-2.5*geo_a, 2.5*geo_a, 100);
Ypoints = linspace(-2.5*geo_a, 2.5*geo_a, 100);
[Hinc, Hscat, Htot] = a.evaluatePECFields(Xpoints, Ypoints, 2);

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
title("Incident Field [A/m]");
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
title("Scattered Field [A/m]");
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
title("Total Field [A/m]");
xlabel("x [m]")
ylabel("y [m]")