% Created by David Newman, PeopleSoft ID 1441532

function RHS = CRS_first_half_step(u,x_interior_points,y_interior_points,LX,LY,delY,bottom_BC,top_BC,left_BC,right_BC)
% Creates right-hand side for each row in the first half step, each
% row involves a tri-diagonal system

% This array will store the right hand side for an entire sweep through the
% array. The right hand sides are stored as columns.

RHS = zeros(x_interior_points,y_interior_points + 2);

tmp1 = zeros(x_interior_points,1);
tmp2 = zeros(x_interior_points,1);
tmp3 = zeros(x_interior_points,1);

% Right hand side for first row
for i = 1:x_interior_points
    tmp1(i) = (1 - 2*LY)*u(1,i) + 2*LY*u(2,i) - 2*LY*delY*bottom_BC(i);
end

% The first and last element of the bottom row RHS has boundary conditions
tmp1(1) = tmp1(1) + LX*left_BC(1);
tmp1(end) = tmp1(end) + LX*right_BC(1);

RHS(:,1) = tmp1;

% Right hand side for the general interior rows
for j = 2:y_interior_points + 1
    for k = 2:x_interior_points + 1
        tmp2(k-1) = LY*u(j-1,k-1) + (1 - 2*LY)*u(j,k-1) + LY*u(j+1,k-1);
    end
    tmp2(1) = tmp2(1) + LX*left_BC(j);
    tmp2(end) = tmp2(end) + LX*right_BC(j);
    RHS(:,j) = tmp2;
    tmp2 = zeros(x_interior_points,1);
end

% Right hand side for last row
for m = 1:x_interior_points
    tmp3(m) = 2*LY*u(end-1,m) + (1 - 2*LY)*u(end,m) + 2*LY*delY*top_BC(m);
end
tmp3(1) = tmp3(1) + LX*left_BC(end);
tmp3(end) = tmp3(end) + LX*right_BC(end);

RHS(:,end) = tmp3;