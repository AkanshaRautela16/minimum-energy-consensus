% Clear workspace and command window
clear; clc; close all;
hold on
% 1. Prompt the user to input the optimal time t (from the Python script)
t_opt = input('Enter the optimal intersection time t (e.g., from Python): ');

if isempty(t_opt) || t_opt <= 0
    error('Time must be a strictly positive number.');
end

disp(['Plotting 3D Reachable Sets at t = ', num2str(t_opt), 's...']);

% 2. Define the 3rd-Order System Matrices
A = [ 0.01,  1.0,  0.0; 
     -1.0,   0.01, 1.0;
      0.0,  -1.0, -1.0];
B = [0; 
     0; 
     1];
n = size(A, 1);

% 3. Compute State Transition Matrix e^(At) and Gramian Wc
M = [A,            B*B'; 
     zeros(n, n), -A'];

E = expm(M * t_opt);
eAt = E(1:n, 1:n);
Wc = E(1:n, n+1:2*n) * E(n+1:2*n, n+1:2*n)';

Wc_inv = inv(Wc);
Wc_half = real(sqrtm(Wc));

% 4. Define Initial Conditions (Columns are the 4 different x0 vectors)
x0_list = [0, 1, 0, 0; 
           0, 0, 1, 0;
           0, 0, 0, 1]; 

% Compute the 4 Ellipsoid Centers
centers = eAt * x0_list;

% 5. Find the Exact Intersection Coordinate using Minimax (fminsearch)
% Objective: Find the point 'x' that minimizes the maximum energy among all 4 ellipsoids
max_energy = @(x) max([ ...
    (x - centers(:,1))' * Wc_inv * (x - centers(:,1)), ...
    (x - centers(:,2))' * Wc_inv * (x - centers(:,2)), ...
    (x - centers(:,3))' * Wc_inv * (x - centers(:,3)), ...
    (x - centers(:,4))' * Wc_inv * (x - centers(:,4)) ]);

% Initial guess is the geometric centroid of all centers
c_mean = mean(centers, 2);
options = optimset('Display', 'off');
intersect_pt = fminsearch(max_energy, c_mean, options);

fprintf('Intersection Point Found: x1=%.4f, x2=%.4f, x3=%.4f\n', ...
    intersect_pt(1), intersect_pt(2), intersect_pt(3));
% 3.5 Calculate and Print Individual Energies
Wc = Wc_half * Wc_half';   % Reconstruct the Gramian
Wc_inv = inv(Wc);          % Inverse Gramian for energy metric
x_star = intersect_pt
final_energies = zeros(1, 4);
fprintf('\n--- Final Energy Distribution at t = %.4f s ---\n', t_opt);
for i = 1:4
    % Distance from the ellipsoid center to the optimal point
    diff_vec = centers(:,i) - x_star; 
    
    % Energy metric: d' * Wc^-1 * d
    final_energies(i) = diff_vec' * Wc_inv * diff_vec; 
    
    % Print the result mapped to its initial condition
    fprintf('Agent %d (Start: [%.0f, %.0f, %.0f]): %.6f\n', ...
        i, x0_list(1,i), x0_list(2,i), x0_list(3,i), final_energies(i));
end

%% 5. Simulate Optimal State Trajectories
disp('Simulating optimal trajectories...');

% 5.1 Precompute constants for the optimal control law
% We need the exact Gramian and state transition at t_opt
M_opt = [A, B*B'; zeros(n, n), -A'];
E_opt = expm(M_opt * t_opt);
eAt_opt = E_opt(1:n, 1:n);

Wc_opt = E_opt(1:n, n+1:2*n) * eAt_opt';
Wc_opt = (Wc_opt + Wc_opt') / 2; % Symmetrize
Wc_inv = inv(Wc_opt);

% 5.2 Compute the constant costate multiplier d_i for each agent
% u_i(t) = B' * expm(A' * (t_opt - t)) * d_i
d = zeros(n, 4);
for i = 1:4
    c_i = eAt_opt * x0_list(:, i);       % Autonomous drift center
    d(:, i) = Wc_inv * (x_star - c_i);    % Costate vector
end

% 5.3 Setup the ODE Simulation
t_span = [0, t_opt];

figure('Color', 'w', 'Position', [200, 200, 800, 650]);
hold on; grid on; view(3);
colors = lines(4);

% Plot the starting positions (Initial Conditions)
for i = 1:4
    plot3(x0_list(1,i), x0_list(2,i), x0_list(3,i), 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end

% 5.4 Simulate and plot each agent
for i = 1:4
    % Define the continuous-time optimal control input u*(t)
    u_opt = @(t) B' * expm(A' * (t_opt - t)) * d(:, i);
    
    % Define the closed-loop ODE: dx/dt = A*x + B*u*(t)
    ode_fun = @(t, x) A*x + B * u_opt(t);
    
    % Integrate the dynamics from t=0 to t_opt
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_sim, x_sim] = ode45(ode_fun, t_span, x0_list(:, i), options);
    
    % Plot the resulting 3D trajectory path
    plot3(x_sim(:,1), x_sim(:,2), x_sim(:,3), 'Color', colors(i,:), ...
        'LineWidth', 2.5, 'DisplayName', ['$x_ ', num2str(i), '(t)$']);
end

% 5.5 Plot the Rendezvous Point
plot3(x_star(1), x_star(2), x_star(3), 'kp', 'MarkerSize', 18, 'DisplayName', 'Optimal Rendezvous $\bar{\mathbf{x}}$');

% Formatting
xlabel('$x_1(t)$', 'Interpreter', 'latex'); 
ylabel('$x_2(t)$', 'Interpreter', 'latex'); 
zlabel('$x_3(t)$', 'Interpreter', 'latex');
%title(['\textbf{Minimum-Energy Trajectories to Rendezvous at } $t = ', num2str(t_opt, '%.4f'), '$s'], ...
 %   'Interpreter', 'latex', 'FontSize', 15);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
axis equal; camlight; lighting gouraud;

% Verify the Maximum
fprintf('----------------------------------------\n');
fprintf('Maximum Energy (Should be exactly 1): %.6f\n\n', max(final_energies));

% 6. Prepare the 3D Plot
figure('Name', sprintf('3D Reachable Sets at t = %.4f', t_opt), 'Color', 'w', 'Position', [100, 100, 900, 800]);
hold on; grid on;
view(3); % Set to 3D view

% Generate a standard 3D unit sphere
[X_sphere, Y_sphere, Z_sphere] = sphere(50);
sphere_pts = [X_sphere(:)'; Y_sphere(:)'; Z_sphere(:)'];

colors = lines(4); % Generate 4 distinct colors
labels = {'Origin (0,0,0)', 'Offset x_1', 'Offset x_2', 'Offset x_3'};

% 7. Transform and Plot each Ellipsoid
for i = 1:4
    c = centers(:, i);
    
    % Affine transformation: x = center + Wc_half * unit_sphere
    ellipsoid_pts = c + Wc_half * sphere_pts;
    
    % Reshape back into 2D grids for surf plotting
    X_ell = reshape(ellipsoid_pts(1,:), size(X_sphere));
    Y_ell = reshape(ellipsoid_pts(2,:), size(Y_sphere));
    Z_ell = reshape(ellipsoid_pts(3,:), size(Z_sphere));
    
    % Plot the 3D surface with 20% transparency (FaceAlpha)
    surf(X_ell, Y_ell, Z_ell, 'FaceColor', colors(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.2, 'DisplayName', labels{i});
    
    % Plot the center point
    plot3(c(1), c(2), c(3), 'x', 'MarkerSize', 10, 'LineWidth', 2, ...
        'Color', colors(i,:), 'HandleVisibility', 'off');
end

% 8. Plot the True Intersection Point in Bold Black
plot3(intersect_pt(1), intersect_pt(2), intersect_pt(3), 'p', ...
    'MarkerSize', 18, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Intersection Point');

% Draw a dashed line from origin to intersection point
plot3([0, intersect_pt(1)], [0, intersect_pt(2)], [0, intersect_pt(3)], ...
    'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% 9. Aesthetics and Lighting
title(sprintf('3D Reachable Ellipsoids Intersecting at t = %.4f s', t_opt), 'FontSize', 14);
xlabel('x_1', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('x_2', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('x_3', 'FontSize', 12, 'FontWeight', 'bold');

axis equal;
camlight; lighting gouraud; % Adds beautiful 3D shading/shadows
legend('Location', 'best', 'FontSize', 11);

% Enable 3D rotation by default
rotate3d on;