% Created by David Newman, PeopleSoft ID 1441532
% Solves the two-dimensional diffusion equation with two Neumann boundary
% conditions using the explicit discretization

clear;
clc;
close all;

% Obtain parameters needed to solve equation.
[x_interior_points,y_interior_points,t_steps,ax,ay,bx,by,T_max,x,y,t,bottom_BC,top_BC,left_BC,right_BC, init, D] = Parameters();

% Lambda x and lambda y appear frequently in the computation.
delT = t(2) - t(1);
delX = x(2) - x(1);
delY = y(2) - y(1);
LX = D*delT/(delX^2);
LY = D*delT/(delY^2);
A = 1 - 2*LX - 2*LY;

if (LX + LY) > 0.5
    fprintf(2,'Stability condition NOT met, stopping execution\n');
    return;
end

% Initialize computational domain. Pad exterior with left and right
% boundary conditions for ease of computation.
u = zeros(y_interior_points+2,x_interior_points + 2,length(t));
for w = 1:length(t)
    u(:,:,w) = [left_BC', zeros(y_interior_points+2,x_interior_points), right_BC'];
end

% Main loop: Loop through each row and determine the solution at the new
% time step. The first and last rows are treated specially due to the
% Neumann boundary condition at the top and bottom of the domain.
for r = 2:length(t)
    u_t = u(:,:,r-1);
    for i = 2:x_interior_points + 1
        u(1,i,r) = A*u_t(1,i) + LX*(u_t(1,i-1) + u_t(1,i+1)) + LY*(2*u_t(2,i) - 2*delY*bottom_BC(i));
        u(end,i,r) = A*u_t(end,i) + LX*(u_t(end,i-1) + u_t(end,i + 1)) + LY*(2*u_t(end-1,i) + 2*delY*top_BC(i));
    end
    for j = 2:y_interior_points + 1
        for k = 2:x_interior_points + 1
            u(j,k,r) = A*u_t(j,k) + LX*(u_t(j,k-1) + u_t(j,k+1)) + LY*(u_t(j-1,k) + u_t(j+1,k));
        end
    end
    fprintf('Time step %d\n',r);
end

% Only keep interior points
u = u(2:end-1,2:end-1,:);

% To verify that the solution has reached steady state
lastChange = u(:,:,end) - u(:,:,end-1);
fprintf('Maximum change in last two time steps = %.15f\n',max(abs(lastChange(:))));

[X,Y] = meshgrid(x(2:end-1),y(2:end-1));
animate(u,X,Y,length(t));
surf(X,Y,u(:,:,end));


