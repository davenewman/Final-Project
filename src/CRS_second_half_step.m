% Created by David Newman, PeopleSoft ID 1441532

function RHS = CRS_second_half_step(u,x_interior_points,y_interior_points,LX,LY,delY,bottom_BC,top_BC,left_BC,right_BC)
% Creates right-hand side for each row in the second half step, each
% column involves a tri-diagonal system

% This array will store the right hand side for an entire sweep through the
% array. The right hand sides are stored as columns.

RHS = zeros(y_interior_points + 2,x_interior_points);

tmp1 = zeros(y_interior_points + 2,1);
tmp2 = zeros(y_interior_points + 2,1);
tmp3 = zeros(y_interior_points + 2,1);

% Right hand side for first column
for i = 1:y_interior_points + 2
    tmp1(i) = (1 - 2*LX)*u(i,1) + LX*u(i,2);
end
tmp1(1) = tmp1(1) - 2*LY*delY*bottom_BC(1);
tmp1(end) = tmp1(end) - 2*LY*delY*top_BC(1);
tmp1 = tmp1 + LX*left_BC';

RHS(:,1) = tmp1;


% Right hand side for general interior columns
for j = 2:x_interior_points - 1
    for k = 1:y_interior_points + 2
        tmp2(k) = LX*u(k,j-1) + (1 - 2*LX)*u(k,j) + LX*u(k,j+1);
    end
    tmp2(1) = tmp2(1) + 2*LY*delY*bottom_BC(j);
    tmp2(end) = tmp2(end) + 2*LY*delY*top_BC(j);
    RHS(:,j) = tmp2;
    tmp2 = zeros(y_interior_points + 2,1);
end

% Right hand side for last column
for m = 1:y_interior_points + 2
    tmp3(m) = LX*u(m,end-1) + (1 - 2*LX)*u(m,end);
end
tmp3(1) = tmp3(1) - 2*LY*delY*bottom_BC(end);
tmp3(end) = tmp3(end) + 2*LY*delY*top_BC(end);
tmp3 = tmp3 + LX*right_BC';
RHS(:,end) = tmp3;
    