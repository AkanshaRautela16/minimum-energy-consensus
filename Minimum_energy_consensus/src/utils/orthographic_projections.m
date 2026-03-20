% Clear workspace and command window
clear; clc; close all;

% 1. Prompt the user to input the optimal time t 
t_opt = input('Enter the optimal intersection time t (e.g., from Python): ');

if isempty(t_opt) || t_opt <= 0
    error('Time must be a strictly positive number.');
end

disp('Computing Projections and 3D Intersection...');

% 2. System Matrices and State Transition
A = [ 0.01,  1.0,  0.0; 
     -1.0,   0.01, 1.0;
      0.0,  -1.0, -1.0];
B = [0; 0; 1];
n = size(A, 1);

M = [A, B*B'; zeros(n, n), -A'];
E = expm(M * t_opt);
eAt = E(1:n, 1:n);
Wc = E(1:n, n+1:2*n) * E(n+1:2*n, n+1:2*n)';
Wc_inv = inv(Wc);

% 3. Extract 2D Projection Shape Matrices
Wc_XY = Wc(1:2, 1:2);       % Projection on x1-x2
Wc_XZ = Wc([1,3], [1,3]);   % Projection on x1-x3
Wc_YZ = Wc(2:3, 2:3);       % Projection on x2-x3

Wc_XY_half = real(sqrtm(Wc_XY));
Wc_XZ_half = real(sqrtm(Wc_XZ));
Wc_YZ_half = real(sqrtm(Wc_YZ));

% 4. Centers and Minimax Intersection
x0_list = [0, 1, 0, 0; 
           0, 0, 1, 0;
           0, 0, 0, 1]; 
centers = eAt * x0_list;

max_energy = @(x) max([ ...
    (x - centers(:,1))' * Wc_inv * (x - centers(:,1)), ...
    (x - centers(:,2))' * Wc_inv * (x - centers(:,2)), ...
    (x - centers(:,3))' * Wc_inv * (x - centers(:,3)), ...
    (x - centers(:,4))' * Wc_inv * (x - centers(:,4)) ]);

c_mean = mean(centers, 2);
options = optimset('Display', 'off');
intersect_pt = fminsearch(max_energy, c_mean, options);

% 5. Setup the Figure for Orthographic Layout
figure('Name', sprintf('Orthographic Projections at t = %.4f', t_opt), 'Color', 'w', 'Position', [50, 50, 1200, 900]);

colors = lines(4);
labels = {'Origin', 'Offset x_1', 'Offset x_2', 'Offset x_3'};

% Define the circle for the 2D projections
theta = linspace(0, 2*pi, 100);
circle = [cos(theta); sin(theta)];

% --- Subplot 1: XY Plane (Top View) ---
subplot(2, 2, 1);
plot_2d_projection(centers, Wc_XY_half, intersect_pt, 1, 2, colors, labels, {'$x_1$', '$x_2$'}, circle);
title('$x_1$ v/s $x_2$','FontSize', 12);

% --- Subplot 2: XZ Plane (Front View) ---
subplot(2, 2, 3);
plot_2d_projection(centers, Wc_XZ_half, intersect_pt, 1, 3, colors, labels, {'$x_1$', '$x_3$'}, circle);
title('$x_1$ v/s $x_3$', 'FontSize', 12);

% --- Subplot 3: YZ Plane (Side View) ---
subplot(2, 2, 4);
plot_2d_projection(centers, Wc_YZ_half, intersect_pt, 2, 3, colors, labels, {'$x_2$', '$x_3$'}, circle);
title('$x_2$ v/s $x_3$','FontSize', 12);

% --- Subplot 4: Full 3D View ---
subplot(2, 2, 2);
hold on; grid on; view(3);
[X_sphere, Y_sphere, Z_sphere] = sphere(40);
sphere_pts = [X_sphere(:)'; Y_sphere(:)'; Z_sphere(:)'];
Wc_3D_half = real(sqrtm(Wc));

for i = 1:4
    c = centers(:, i);
    ellipsoid_pts = c + Wc_3D_half * sphere_pts;
    X_ell = reshape(ellipsoid_pts(1,:), size(X_sphere));
    Y_ell = reshape(ellipsoid_pts(2,:), size(Y_sphere));
    Z_ell = reshape(ellipsoid_pts(3,:), size(Z_sphere));
    surf(X_ell, Y_ell, Z_ell, 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');
end
plot3(intersect_pt(1), intersect_pt(2), intersect_pt(3), 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Intersection');

title('3D Perspective', 'FontSize', 12);
xlabel('$x_1$', 'FontWeight', 'bold'); ylabel('$x_2$', 'FontWeight', 'bold'); zlabel('$x_3$', 'FontWeight', 'bold');
axis equal; camlight; lighting gouraud; rotate3d on;

% Add a master legend to the 3D plot using dummy lines for neatness
h1 = plot(nan, nan, 'Color', colors(1,:), 'LineWidth', 2);
h2 = plot(nan, nan, 'Color', colors(2,:), 'LineWidth', 2);
h3 = plot(nan, nan, 'Color', colors(3,:), 'LineWidth', 2);
h4 = plot(nan, nan, 'Color', colors(4,:), 'LineWidth', 2);
hp = plot(nan, nan, 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
legend([h1, h2, h3, h4, hp], [labels, {'Intersection Point'}], 'Location', 'bestoutside');

% Master Title
sgtitle(sprintf('Orthographic Projections of Reachable Sets at t = %.4f s', t_opt), 'FontSize', 16, 'FontWeight', 'bold');


% =========================================================================
% LOCAL FUNCTIONS (Must be at the very bottom of the MATLAB script!)
% =========================================================================
function plot_2d_projection(centers, Wc_half, intersect_pt, idx1, idx2, colors, labels, axis_labels, circle_pts)
    hold on; grid on;
    for i = 1:4
        c_2d = centers([idx1, idx2], i);
        ellipse_pts = c_2d + Wc_half * circle_pts;
        plot(ellipse_pts(1,:), ellipse_pts(2,:), 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', labels{i});
        plot(c_2d(1), c_2d(2), 'x', 'MarkerSize', 8, 'Color', colors(i,:), 'HandleVisibility', 'off');
    end
    % Plot Intersection Point
    plot(intersect_pt(idx1), intersect_pt(idx2), 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Intersection Point');
    xlabel(axis_labels{1}, 'FontWeight', 'bold');
    ylabel(axis_labels{2}, 'FontWeight', 'bold');
    axis equal;
end