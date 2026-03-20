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
        'LineWidth', 2.5, 'DisplayName', ['Agent ', num2str(i)]);
end

% 5.5 Plot the Rendezvous Point
plot3(x_star(1), x_star(2), x_star(3), 'kp', 'MarkerSize', 18, ...
    'MarkerFaceColor', 'y', 'DisplayName', 'Optimal Rendezvous $x^*$');

% Formatting
xlabel('$x_1$', 'Interpreter', 'latex'); 
ylabel('$x_2$', 'Interpreter', 'latex'); 
zlabel('$x_3$', 'Interpreter', 'latex');
title(['\textbf{Minimum-Energy Trajectories to Rendezvous at } $t = ', num2str(t_opt, '%.4f'), '$s'], ...
    'Interpreter', 'latex', 'FontSize', 15);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
axis equal; camlight; lighting gouraud;