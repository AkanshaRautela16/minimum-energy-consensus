% Clear workspace and command window
clear; clc; close all;

% 1. Prompt the user to input the time t (from the SageMath output)
t = input('Enter the optimal time t (e.g., from SageMath): ');

if isempty(t) || t <= 0
    error('Time must be a strictly positive number.');
end

% 2. Define the NEW System Matrices
A = [0.01, 1; 
    -1,    0.01];
B = [0; 
     1];

% 3. Compute State Transition Matrix e^(At)
eAt = expm(A * t);

% 4. Compute Finite-Time Gramian Wc(t) using the Matrix Exponential Trick
n = size(A, 1);
M = [A,            B*B'; 
     zeros(n, n), -A'];

% Compute exponential of the block matrix
E = expm(M * t);

% Extract the sub-blocks to get the Gramian
E12 = E(1:n, n+1:2*n);
E22 = E(n+1:2*n, n+1:2*n);
Wc = E12 * E22';

% 5. Define Initial Conditions
x0_list = [0, 1, 2; 
           0, 0, 0]; 

labels = {'Origin (0,0)', 'Offset x_1 (1,0)', 'Offset x_2 (0,1)'};
colors = lines(3); % Generate 3 distinct colors

% 6. Prepare the Plot
figure('Name', sprintf('Reachable Sets at t = %.4f', t), 'Color', 'w');
hold on; grid on;
xlabel('x_1', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('x_2', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Reachable Ellipsoids at time t = %.4f (E_{max} = 1)', t), 'FontSize', 14);

% Generate 200 points around a standard unit circle
theta = linspace(0, 2*pi, 2000);
unit_circle = [cos(theta); sin(theta)];

% Matrix square root of Wc maps the unit circle exactly onto the boundary
Wc_half = sqrtm(Wc);

% 7. Plot each ellipsoid
for i = 1:size(x0_list, 2)
    x0 = x0_list(:, i);
    
    % Compute the center of the ellipsoid (zero-input response: e^At * x0)
    center = eAt * x0;
    
    % Compute the actual coordinates of the ellipsoid boundary
    ellipsoid_pts = center + Wc_half * unit_circle;
    
    % Plot the boundary curve
    plot(ellipsoid_pts(1,:), ellipsoid_pts(2,:), 'LineWidth', 2, ...
        'Color', colors(i,:), 'DisplayName', labels{i});
    
    % Plot the center dot
    plot(center(1), center(2), '.', 'MarkerSize', 25, 'Color', colors(i,:), ...
        'HandleVisibility', 'off'); 
end

% Ensure the X and Y axes are scaled equally so ellipses aren't distorted
axis equal;
legend('Location', 'best', 'FontSize', 11);