
% Evolution of states with t_f%
clc; clear all; close all;
x01 = 1;  
x02 = 1;
beta = 2;

figure;
hold on;
grid on;

for tf = [5,10,15,25]
    f = @(x1,x2) ([x1 - x01 - tf*x02; x2 - x02]' * ...
                  ([tf^3/3 tf^2/2; tf^2/2 tf])^-1 * ...
                  [x1 - x01 - tf*x02; x2 - x02]) - beta;

    % Plot the implicit function
    fimplicit(f, [-50 70 -10 12], 'DisplayName', ['t_f = ', num2str(tf)]);
end

legend('show');
xlabel('x_1'); ylabel('x_2');


fimplicit(@f,[-50, 70, -10, 12])

