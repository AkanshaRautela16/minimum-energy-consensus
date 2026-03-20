%% 6. Plot Optimal Control Input Profiles
disp('Calculating optimal control inputs...');

% 6.1 Define time vector for plotting
N_pts = 200;
t_vec = linspace(0, t_opt, N_pts);
u_profiles = zeros(4, N_pts); % 4 agents, N_pts time steps

% 6.2 Compute u(t) for each agent across the time horizon
for i = 1:4
    for k = 1:N_pts
        t = t_vec(k);
        % Evaluate the analytical optimal control law
        u_profiles(i, k) = B' * expm(A' * (t_opt - t)) * d(:, i);
    end
end

% 6.3 Visualization
figure('Color', 'w', 'Position', [250, 250, 800, 400]);
hold on; grid on;

% Re-use the same color palette so Agent 1 matches the 3D plot
colors = lines(4); 

for i = 1:4
    plot(t_vec, u_profiles(i, :), 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', ['$u^*_ ', num2str(i), '(t)$']);
end

% 6.4 Formatting
xlabel('Time $t$ (seconds)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Control Inputs', 'Interpreter', 'latex', 'FontSize', 12);
%title('\textbf{Optimal Control Input Profiles } $u^*(t)$', 'Interpreter', 'latex', 'FontSize', 15);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
xlim([0, t_opt]);

% Add a subtle zero-line to show positive/negative actuation
yline(0, 'k--', 'HandleVisibility', 'off'); 

% 6.5 Numerical Verification of Energy Constraint
fprintf('--- Verification: Integral of u^2(t) dt ---\n');
for i = 1:4
    % Numerically integrate the square of the input signal using trapz
    energy_used = trapz(t_vec, u_profiles(i, :).^2);
    fprintf('Agent %d Total Energy Consumed: %.6f\n', i, energy_used);
end
fprintf('----------------------------------------\n\n');