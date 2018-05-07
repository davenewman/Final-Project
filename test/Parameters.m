function [x_interior_points,y_interior_points,t_steps,ax,ay,bx,by,T_max,x,y,t,bottom_BC,top_BC,left_BC,right_BC, init, D] = Parameters()
% This file holds the parameters to be loaded into the main script. Used
% for modularity of code.

% Number of points to discretize the domain
x_interior_points = 3;
y_interior_points = 3;
t_steps = 10000;

% Diffusion coefficient (1 for this problem)
D = 1;

% Domain
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;
T_max = 10;

% Set up vectors for boundaries
x = linspace(ax,bx,x_interior_points+2);
y = linspace(ay,by,y_interior_points+2);
t = linspace(0,T_max,t_steps);

% Bottom BC (Von Neumann)
bottom_BC = zeros(1,length(x));

% Top BC (Von Neumann)
top_BC = zeros(1,length(x));

% Left BC (Dirichlet)
left_BC = (y.^2).*sin(y/4);


% Right BC (Dirichlet)
right_BC = cos(pi*y).*cosh(2*pi-y);


% Initial condition
init = zeros(y_interior_points + 2,x_interior_points);