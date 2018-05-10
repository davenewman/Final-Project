% Created by David Newman, PeopleSoft ID 1441532
% Solves the two-dimensional diffusion equation with two Neumann boundary
% conditions using the ADI (Alternating Direction Implicit) discretization

clear;
clc;
close all;

delta_x = [];
L1h =[];
L2h =[];
Lih =[];

for zz = [10 20 40 80 160 320];

% Number of points to discretize the domain
x_interior_points = zz;
y_interior_points = zz;
t_steps = 1000;

tau = -2;
k = 1;

% Diffusion coefficient (1 for this problem)
D = -1/(2*tau*k^2);

% Domain
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;
T_max = 2;

% Set up vectors for boundaries
x = linspace(ax,bx,x_interior_points+2);
y = linspace(ay,by,y_interior_points+2);
t = linspace(0,T_max,t_steps);

% Bottom BC (Von Neumann)
bottom_BC = zeros(1,length(x));

% Top BC (Von Neumann)
top_BC = zeros(1,length(x));

% Left BC (Dirichlet)
left_BC = zeros(1,length(y));

% Right BC (Dirichlet)
right_BC = zeros(1,length(y));

% Initial condition
init =0*ones(y_interior_points + 2,x_interior_points);
[P,Q] = meshgrid(x(2:end-1),y);
init = sin(P).*cos(Q);

% Obtain parameters needed to solve equation

% Lambda x and lambda y appear frequently in the computation
delT = t(2) - t(1);
delX = x(2) - x(1);
delY = y(2) - y(1);
LX = D*delT/(2*delX^2);
LY = D*delT/(2*delY^2);

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

% Initialize computational domain
u = init;



% Main loop: create right-hand side for first half step, solve the
% tri-diagonal system for each row, create right-hand side for second half
% step, solve the tri-diagonal system for each row.
for r = 2:length(t)
    RHS = CRS_first_half_step(u(:,:,r-1),x_interior_points,y_interior_points,LX,LY,delY,bottom_BC,top_BC,left_BC,right_BC);
    for s = 1:y_interior_points + 2
        u(s,:,r) = SolveTriDiag(a_first_half_step,b_first_half_step,c_first_half_step,RHS(:,s));
    end
    RHS = CRS_second_half_step(u(:,:,r),x_interior_points,y_interior_points,LX,LY,delY,bottom_BC,top_BC,left_BC,right_BC);
    for w = 1:x_interior_points
        u(:,w,r) = SolveTriDiag(a_second_half_step,b_second_half_step,c_second_half_step,RHS(:,w));
    end 
%     fprintf('Time step %d\n',r);
    
end

% Only keep interior points
u = u(2:end-1,:,:);
change = zeros(1,length(t)-1);
% for z = 2:length(t)
%     utemp1 = u(:,:,z-1);
%     utemp2 = u(:,:,z);
%     utemp1 = utemp1(:);
%     utemp2 = utemp2(:);
%     change(z) = max( abs(utemp2 - utemp1));
% end
% semilogy(t,change,'k'),xlabel('Time'),ylabel('Maximum Change'),title('Maximum change in solution between time steps')

[X,Y] = meshgrid(x(2:end-1),y(2:end-1));
% animate_and_save(u,X,Y,length(t));
% surf(X,Y,u(:,:,end/2));
% xlabel('x'),ylabel('y'),title('Solution for n = 80 at T = 2, ADI scheme');

u_exact = exp(-T_max/2).*sin(X).*cos(Y);


erfunc = u(:,:,end);

delta_x = [delta_x delX];
L1h =[L1h  (1/zz^2)*sum( abs( u_exact(:) - erfunc(:) ) )];
L2h =[L2h  (1/zz)*norm( abs( u_exact(:) - erfunc(:) ) )];
Lih =[Lih  max( abs( u_exact(:) - erfunc(:)) )];

fprintf('Iteration\n');

end

figure;
loglog(delta_x,L1h,'ks');
xlabel('Spatial Step (delta x)');
ylabel('Error');
title('L1 error comparison with manufactured solution');
a1 = polyfit(log(delta_x),log(L1h),1)

figure;
loglog(delta_x,L2h,'ks');
xlabel('Spatial Step (delta x)');
ylabel('Error');
title('L2 error comparison with manufactured solution');
a2 = polyfit(log(delta_x),log(L2h),1)

figure;
loglog(delta_x,Lih,'ks');
xlabel('Spatial Step (delta x)');
ylabel('Error');
title('L-Infinity error comparison with manufactured solution');
a3 = polyfit(log(delta_x),log(Lih),1)

