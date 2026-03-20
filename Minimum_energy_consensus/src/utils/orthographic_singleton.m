%% Corrected Minimum Time Singleton Intersection (MATLAB)
clear; clc; close all;

% FIX: Corrected Global LaTeX Interpreters
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); % Note the 'Axes' prefix
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% 1. System Definition
A = [ 0.01,  1.0,  0.0; 
     -1.0,   0.01, 1.0;
      0.0,  -1.0, -1.0];
B = [0; 0; 1];
n = size(A, 1);

% Initial Conditions (Origin + 3 basis vectors)
x0_list = [0, 1, 0, 0; 
           0, 0, 1, 0; 
           0, 0, 0, 1];

% 2. Solve for Optimal Time t
disp('Finding optimal intersection time...');
% Search interval [0.5, 5]
try
    t_opt = fzero(@(t) compute_minimax_energy(t, A, B, x0_list) - 1, [0.5, 5]);
catch
    error('Could not find a root in the interval. Try expanding the search range.');
end

% 3. Compute Final State for plotting
[~, x_star, centers, Wc_half] = compute_minimax_energy(t_opt, A, B, x0_list);

% 4. Visualization
figure('Color', 'w', 'Position', [100, 100, 1100, 850]);
colors = lines(4);

% --- Subplot 1: x1 vs x2 ---
subplot(2,2,1); hold on; grid on;
plot_proj(centers, Wc_half, x_star, 1, 2, colors, '$x_1$', '$x_2$', 'Projection: $x_1$ vs $x_2$');

% --- Subplot 2: 3D View ---
subplot(2,2,2); hold on; grid on; view(3);
[Xs, Ys, Zs] = sphere(40);
sphere_pts = [Xs(:)'; Ys(:)'; Zs(:)'];
for i = 1:4
    pts = centers(:,i) + Wc_half * sphere_pts;
    surf(reshape(pts(1,:), size(Xs)), reshape(pts(2,:), size(Ys)), reshape(pts(3,:), size(Zs)), ...
        'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');
end
plot3(x_star(1), x_star(2), x_star(3), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', '$x^*$');
xlabel('$x_1$'); ylabel('$x_2$'); zlabel('$x_3$');
title('\textbf{3D State Space Reachable Sets}');
axis equal; camlight; lighting gouraud;
legend('Location', 'northeast');

% --- Subplot 3: x1 vs x3 ---
subplot(2,2,3); hold on; grid on;
plot_proj(centers, Wc_half, x_star, 1, 3, colors, '$x_1$', '$x_3$', 'Projection: $x_1$ vs $x_3$');

% --- Subplot 4: x2 vs x3 ---
subplot(2,2,4); hold on; grid on;
plot_proj(centers, Wc_half, x_star, 2, 3, colors, '$x_2$', '$x_3$', 'Projection: $x_2$ vs $x_3$');

sgtitle(['\textbf{Singleton Intersection at } $t = ', num2str(t_opt, '%.4f'), '$s'], 'FontSize', 16);

%% --- Helper Functions ---

function [energy, x_opt, centers, Wc_half] = compute_minimax_energy(t, A, B, x0_list)
    n = size(A, 1);
    % Van Loan's Method for Gramian
    M = [A, B*B'; zeros(n, n), -A'];
    E = expm(M * t);
    eAt = E(1:n, 1:n);
    % Extract Wc and ensure strict symmetry
    Wc = E(1:n, n+1:2*n) * eAt'; 
    Wc = (Wc + Wc') / 2; 
    
    Wc_inv = inv(Wc);
    Wc_half = real(sqrtm(Wc));
    centers = eAt * x0_list;
    
    % Objective function
    obj = @(x) max_calc(x', centers, Wc_inv);
    
    options = optimset('Display', 'off', 'TolX', 1e-7);
    [x_opt, energy] = fminsearch(obj, mean(centers, 2), options);
    x_opt = x_opt';
end

function m = max_calc(x, centers, Wc_inv)
    diffs = centers - x; 
    energies = zeros(1, 4);
    for i = 1:4
        energies(i) = diffs(:,i)' * Wc_inv * diffs(:,i);
    end
    m = max(energies);
end

function plot_proj(centers, Wc_half, x_star, idx1, idx2, colors, lab1, lab2, titl)
    theta = linspace(0, 2*pi, 200);
    circle = [cos(theta); sin(theta)];
    Wc_full = Wc_half * Wc_half';
    Wc_proj = Wc_full([idx1, idx2], [idx1, idx2]);
    Wp_half = real(sqrtm(Wc_proj));
    for i = 1:4
        pts = centers([idx1, idx2], i) + Wp_half * circle;
        plot(pts(1,:), pts(2,:), 'LineWidth', 1.8, 'Color', colors(i,:), 'HandleVisibility', 'off');
    end
    plot(x_star(idx1), x_star(idx2), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
    xlabel(lab1); ylabel(lab2); title(titl);
    axis equal; grid on;
end