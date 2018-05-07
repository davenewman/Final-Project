% Created by David Newman, PeopleSoft ID 1441532

clear;
clc;
close all;

% Obtain parameters needed to solve equation. For now we are only dealing
% with the ADI method however there may be some modularity added later to
% deal with the explicit method, which should be MUCH easier to implement.

[x_interior_points,y_interior_points,t_steps,ax,ay,bx,by,T_max,x,y,t,bottom_BC,top_BC,left_BC,right_BC, init, D] = Parameters();

% Lambda x and lambda y appear frequently in the computation
delT = t(2) - t(1);
delX = x(2) - x(1);
delY = y(2) - y(1);
LX = D*delT/(2*delX^2);
LY = D*delT/(2*delY^2);

% n is the total number of interior points the "+2" comes from the extra
% unknowns due to the Neumann boundary conditions on top and bottom
n = x_interior_points*(y_interior_points + 2);

% Build the diagonals of the main array at each half step. "a" is the main
% diagonal, "b" is the lower diagonal, and "c" is the upper diagonal. For
% the first half step, a tri-diagonal system is solved for each row of the 
% domain. For the second half step, a tri-diagonal system is solved for 
% each column of the domain.
a_first_half_step = (1 + 2*LX)*ones(1,x_interior_points);
b_first_half_step = -LX*ones(1,x_interior_points - 1);
c_first_half_step = b_first_half_step;

a_second_half_step = (1 + 2*LY)*ones(1,y_interior_points + 2);
b_second_half_step = [-LY*ones(1,y_interior_points) -2*LY];
c_second_half_step = [-2*LY -LY*ones(1,y_interior_points)];

u = init;

% TODO: Main loop:
        % /DONE * Create right side for initial half step (implicit in x) 
        % * Solve tri-diagonal system for first half step
for r = 2:length(t)
    RHS = CRS_first_half_step(u(:,:,r-1),x_interior_points,y_interior_points,LX,LY,delY,bottom_BC,top_BC,left_BC,right_BC);
    for s = 1:y_interior_points + 2
        u(s,:,r) = SolveTriDiag(a_first_half_step,b_first_half_step,c_first_half_step,RHS(:,s));
    end
    RHS = CRS_second_half_step(u(:,:,r),x_interior_points,y_interior_points,LX,LY,delY,bottom_BC,top_BC,left_BC,right_BC);
    for w = 1:x_interior_points
        u(:,w,r) = SolveTriDiag(a_second_half_step,b_second_half_step,c_second_half_step,RHS(:,w));
    end 
end

u = u(2:end-1,:,:);

lastError = u(:,:,end) - u(:,:,end-1)

surf(u(:,:,end))
xlabel('x');
ylabel('y')



