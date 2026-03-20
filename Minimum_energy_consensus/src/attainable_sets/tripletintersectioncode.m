clc; clear all; close all;
[0.04, 0.1; 0.39, 1.05; 0.3, -0.525; 0.5, -0.525]
% Parameters
x0  = [0.04,  0.10];
x1  = [0.39,  1.05];
x2  = [-0.9, -0.72];
beta = 0.7;
tf   = 6.23;

% choose colors 
C0 = [0.8500 0.3250 0.0980];
C1 = [0.0000 0.4470 0.7410];
C2 = [0.9290 0.6940 0.1250];

% Shared matrix and quadratic form (stable solve, no inv)
Q = [tf^3/3, tf^2/2; tf^2/2, tf];                               
quad = @(dx) sum((Q \ dx).*dx,1) - beta;                           
%quad = @(dx) dx'*Q^-1*dx - beta;
% Implicit functions: f(x,y)=0
f0 = @(x,y) quad([x - x0(1) - tf*x0(2); y - x0(2)]);               
f1 = @(x,y) quad([x - x1(1) - tf*x1(2); y - x1(2)]);               
f2 = @(x,y) quad([x - x2(1) - tf*x2(2); y - x2(2)]);               

% Figure and axes
figure; ax = axes; hold(ax,'on'); axis(ax,'equal'); grid(ax,'on'); 

box  = [-5 8 -5 5];
dens = 150;                                                        

% Curves with explicit colors
h0 = fimplicit(ax, f0, box, 'MeshDensity', dens, 'Color', C0, 'DisplayName','set x0'); 
h1 = fimplicit(ax, f1, box, 'MeshDensity', dens, 'Color', C1, 'DisplayName','set x1');
h2 = fimplicit(ax, f2, box, 'MeshDensity', dens, 'Color', C2, 'DisplayName','set x2');

% Initial-condition markers with matching colors
plot(x0(1), x0(2), '*', 'Color', C0, 'LineWidth', 1.2, 'DisplayName','x0');            
plot(x1(1), x1(2), '*', 'Color', C1, 'LineWidth', 1.2, 'DisplayName','x1');            
plot(x2(1), x2(2), '*', 'Color', C2, 'LineWidth', 1.2, 'DisplayName','x2');            

% Labels and legend
xlabel('$x(t)$', 'FontSize', 14, 'FontWeight','bold','Interpreter','latex');
ylabel('$\dot{x}(t)$', 'FontSize', 14, 'FontWeight','bold','Interpreter','latex');                      