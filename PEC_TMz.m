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

% assign an incidence
a = a.constructIncidence(0, 1);

% solve the method-of-moments system
a = a.solveI();

% plot the mesh
figure();
a.plotMesh();
set(gcf, 'Color', 'w');
set(gca, 'FontName', 'times', 'FontSize', 18, 'TickLength', [0,0])
box on; grid on;
axis equal tight;
title("Mesh")
xlabel("x [m]")
ylabel("y [m]")

% plot the current
figure();

set(gcf, 'Color', 'w');
set(gca, 'FontName', 'times', 'FontSize', 18);
axis equal tight;
title("Surface Current");