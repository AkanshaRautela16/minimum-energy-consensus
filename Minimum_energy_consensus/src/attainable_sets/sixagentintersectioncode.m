clc; clear all; close all;

% Parameters
xc = -0.5 ;  
yc = 0;
IC = [3,-1.1; 
      -3, 1.1;  
      3.1, -3;
      -3.1 ,3; 
       6.5, 6.5 ;
       -7.5,-6.5];
x0 = IC(1,:); 
x1 = IC(2,:);  
x2 = IC(3,:);
x3 = IC(4,:); 
x4 = IC(5,:);
x5 = IC(6,:);
beta = 30;
tf   = 8.158576;

% choose colors 
C0 = [0.8500 0.3250 0.0980];
C1 = [0 0.4470 0.7410];
C2 = [0.9290 0.6940 0.1250];
C3 = [0.4660 0.6740 0.1880];
C4 = [0.3010 0.7450 0.9330];
C5 = [0.6350 0.0780 0.1840];

% Shared matrix and quadratic form (stable solve, no inv)
Q = [tf^3/3, tf^2/2; tf^2/2, tf];                               
quad = @(dx) sum((Q \ dx).*dx,1) - beta;                           
%quad = @(dx) dx'*Q^-1*dx - beta;
% Implicit functions: f(x,y)=0
f0 = @(x,y) quad([x - x0(1) - tf*x0(2); y - x0(2)]);               
f1 = @(x,y) quad([x - x1(1) - tf*x1(2); y - x1(2)]);               
f2 = @(x,y) quad([x - x2(1) - tf*x2(2); y - x2(2)]); 
f3 = @(x,y) quad([x - x3(1) - tf*x3(2); y - x3(2)]); 
f4 = @(x,y) quad([x - x4(1) - tf*x4(2); y - x4(2)]); 
f5 = @(x,y) quad([x - x5(1) - tf*x5(2); y - x5(2)]); 
axis  = [-52 52 -5 5];
% Figure and axes
figure; ax = axes; hold(ax,'on'); grid(ax,'on'); 

box  = [-135 135 -25 25];
dens = 150;                                                        

% Curves with explicit colors
h0 = fimplicit(ax, f0, box, 'MeshDensity', dens, 'Color', C0, 'LineWidth', 2, 'DisplayName','set x0'); 
h1 = fimplicit(ax, f1, box, 'MeshDensity', dens, 'Color', C1, 'LineWidth', 2, 'DisplayName','set x1');
h2 = fimplicit(ax, f2, box, 'MeshDensity', dens, 'Color', C2, 'LineWidth', 2, 'DisplayName','set x2');
h3 = fimplicit(ax, f3, box, 'MeshDensity', dens, 'Color', C3, 'LineWidth', 2, 'DisplayName','set x3');
h4 = fimplicit(ax, f4, box, 'MeshDensity', dens, 'Color', C4, 'LineWidth', 2, 'DisplayName','set x4');
h5 = fimplicit(ax, f5, box, 'MeshDensity', dens, 'Color', C5, 'LineWidth', 2, 'DisplayName','set x5');

% Initial-condition markers with matching colors
plot(x0(1), x0(2), '+', 'Color', C0, 'LineWidth', 1.5, 'DisplayName','x0');            
plot(x1(1), x1(2), '+', 'Color', C1, 'LineWidth', 1.5, 'DisplayName','x1');            
plot(x2(1), x2(2), '+', 'Color', C2, 'LineWidth', 1.5, 'DisplayName','x2');  
plot(x3(1), x3(2), '+', 'Color', C3, 'LineWidth', 1.5, 'DisplayName','x3'); 
plot(x4(1), x4(2), '+', 'Color', C4, 'LineWidth', 1.5, 'DisplayName','x4'); 
plot(x5(1), x5(2), '+', 'Color', C5, 'LineWidth', 1.5, 'DisplayName','x5'); 
plot(xc,yc,'p', 'MarkerSize', 5, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'y', 'LineWidth', 2)
% Labels and legend
xlabel('$x(t)$', 'FontSize', 14, 'FontWeight','bold','Interpreter','latex');
ylabel('$\dot{x}(t)$', 'FontSize', 14, 'FontWeight','bold','Interpreter','latex');                      