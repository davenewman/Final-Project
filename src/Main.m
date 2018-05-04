% Created by David Newman, PeopleSoft ID 1441532

clear;
clc;
close all;

% Obtain parameters needed to solve equation. For now we are only dealing
% with the ADI method however there may be some modularity added later to
% deal with the explicit method, which should be MUCH easier to implement.

[x_interior_points,y_interior_points,t_steps,ax,ay,bx,by,T_max,x,y,t,bottom_BC,top_BC,left_BC,right_BC] = Parameters();

% TODO: Handle Neumann BCs and turn them into Dirichlet with second order
% formula

% TODO: Main loop:
        % * Create right side for initial half step (implicit in x)
        % * Solve tri-diagonal system for first half step
        % * Solve tri-diagonal system for second half step (same
        % coefficients, implicit in y)
        % * Store solution in 2D array? Initialize before loop start
        % * Advance time step
        
% TODO: Reshape 2D array into 3D array, each page representing a different
% time step

% TODO: Solve explicit method, hopefully with this same script. Ensure
% check on stability condition is carried out
