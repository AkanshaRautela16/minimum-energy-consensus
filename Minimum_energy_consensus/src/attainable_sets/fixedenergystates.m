clc; clear; close all;

x01 = 1;  
x02 = 1;  
beta = 2;

figure; hold on; grid on;

xlim([-50 70])
ylim([-10 15])

for tf = [5,10,15,25]

    A = [tf^3/3 tf^2/2; tf^2/2 tf];
    Ainv = inv(A);

    f = @(x1,x2) ([x1 - x01 - tf*x02; x2 - x02]' * Ainv * ...
                  [x1 - x01 - tf*x02; x2 - x02]) - beta;

    fimplicit(f, [-50 70 -10 15], ...
        'DisplayName', ['t_f = ', num2str(tf)]);
end
plot(x01,x02,'+', 'MarkerSize', 8, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'y', 'LineWidth', 2)
legend show
xlabel('$x(t)$', 'FontSize', 16, 'fontweight','bold','Interpreter','latex');
ylabel('$\dot{x}(t)$', 'Fontsize', 16, 'fontweight','bold','Interpreter','latex');