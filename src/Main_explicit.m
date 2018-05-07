% Created by David Newman, PeopleSoft ID 1441532
% Solves the two-dimensional diffusion equation with two Neumann boundary
% conditions using an explicit discretization

clear;
close all;
clc;

% Obtain parameters needed to solve equation.
[x_interior_points,y_interior_points,t_steps,ax,ay,bx,by,T_max,x,y,t,bottom_BC,top_BC,left_BC,right_BC, init, D] = Parameters();

delT = t(2) - t(1);
delX = x(2) - x(1);
delY = y(2) - y(1);
LX = D*delT/(delX^2);
LY = D*delT/(delY^2);
A = 1 - 2*LX - 2*LY;

if (LX + LY) > 0.25
    fprintf(2,'Stability condition NOT met\n');
end

% No need for expensive function calls here, this main loop should be clear
% enough. We will iterate through the solution given by the previous time
% step, treating the first and last rows separately as they will have the
% Neumann conditions and their calculation will be different. The
% initialization this time will include the left and right boundaries.

init = [left_BC', zeros(y_interior_points+2,x_interior_points), right_BC'];
u = zeros(y_interior_points + 2, x_interior_points + 2, length(t));
u(:,:,1) = init;

for a = 2:length(t)
    u_t = u(:,:,a-1);
    tmp = zeros(size(init));
    % First row and last row
    for i = 2:x_interior_points + 1
        %tmp(1,i) = % something having to do with bottom boundary
        %tmp(end,i) = % something having to do with top boundary
    end
    
    % Loop through main part of array
    for j = 2:y_interior_points + 1
        for k = 2:x_interior_points + 1
            tmp(j,k) = A*u(j,k) + LX*(u(j,k-1) + u(j,k+1)) + LY*(u(j-1,k) + u(j+1,k));
        end
    end
end
    
            