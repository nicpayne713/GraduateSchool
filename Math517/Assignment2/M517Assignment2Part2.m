% Math 517 Assignment 2 Part 1
clear all
close all
% Construct sparse coefficient matrix A and appropriate right-hand side
% vector F for the descretized version utilizing the 9-point Laplacian on a
% uniform mesh of Poisson's equation in 2D 
% First I will set up the grid
a = 0;
b = 1;
% The grid will be based on the desired step size for the problem
% For the purposes of doing the convergence study, h is a vector and the
% script loops through the problem for how ever many h's there are to try. 
h = [1/2,1/4,1/16];
Error = zeros(1,length(h));
for i = 1:length(h)
    N = (1/h(i))-1;
    % This creates the proper A matrix for (-)Laplacian(u) = f. No sign
    % changes are necessary
    e = ones(N,1);
    S = spdiags([e 10*e e], [-1 0 1], N, N);
    I = spdiags([-1/2*e e -1/2*e], [-1 0 1], N, N);
    A = kron(I,S)+kron(S,I);
    % Set up the interior and boundary points for the grid
    x = linspace(a,b,N+2);
    y = linspace(a,b,N+2);
    [X,Y] = meshgrid(x,y); % 2D cartesian grid of the x,y values
    X = X';Y = Y'; % take transpose so that X(i,j),Y(i,j) is the (i,j)th point.
    
    % Now to set up notating the interior points in x and y
    IntPoints = 2:N+1;
    Xinterior = X(IntPoints,IntPoints);
    Yinterior = Y(IntPoints,IntPoints);
    
    % The test problem
    f = @(x,y) 1.25.*exp(x+.5.*y);
    
    % Uhat is the exact solution sampled on the whole mesh
    Uhat = exp(X+0.5.*Y);
    usol = Uhat;
    
    % f needs to be evaluated at all the interior points, so that's done here

    RHS = 6*h(i)^2*f(Xinterior,Yinterior);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Finding the Laplacian(f) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Lf is h^2/12 * 5pointlaplacian(f). The h^2s cancel out
    Lf = zeros(N,N);
    for j = 1:N
        for k = 1:N
            Lf(j,k) = (f(j-1,k)+f(j+1,k)+f(j,k+1)+ f(j,k-1)-4*f(j,k));
        end
    end
    
    % Change the boundary conditions on the 5 point laplacian of f
%     Lf(:,1) = Lf(:,1) + usol(IntPoints,1)/(h(i)^2);
%     Lf(:,N) = Lf(:,N) + usol(IntPoints,N+2)/(h(i)^2);
%     Lf(1,:) = Lf(1,:) + usol(1,IntPoints)/(h(i)^2);
%     Lf(N,:) = Lf(N,:) + usol(N+2,IntPoints)/(h(i)^2);
    
    % Scale Lf to account for the various coefficients
    Lf = (h(i)^2/2)*Lf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Now we set the boundary conditions, since we have a known solution around
    % the boundary we just need a full array and we'll use the boundary
    % values only

    RHS(:,1) = RHS(:,1) - usol(IntPoints,1) - usol(IntPoints-1,1) - usol(IntPoints+1,1) ;
    RHS(:,N) = RHS(:,N) - usol(IntPoints,N+2)- usol(IntPoints-1,N+2) - usol(IntPoints+1,N+2) ;
    RHS(1,:) = RHS(1,:) - usol(1,IntPoints)- usol(1,IntPoints-1) - usol(1,IntPoints+1) ;
    RHS(N,:) = RHS(N,:) - usol(N+2,IntPoints)- usol(N+2,IntPoints-1) - usol(N+2,IntPoints+1);
    
    % The corner pieces of the interior mesh get the boundary conditions at the
    % corners subtracted twice, so here we add 1 back in.
    RHS(1,1) = RHS(1,1) + usol(1,1);
    RHS(1,N) = RHS(1,N) + usol(1,N+2);
    RHS(N,N) = RHS(N,N) + usol(N+2,N+2);
    RHS(N,1) = RHS(N,1) + usol(N+2,1);
    
    % We need to add the laplacian of the right-hand side function
    RHS = (RHS + Lf);
    
    % Now I can take the grid created with the meshgrid function and stretch it
    % out into an array for the RHS of the system
    F = -reshape(RHS,N*N,1);
    
    % Now to solve the system using the backslash
    U = A\F;
    
    % The u vector gets put into the usol grid as the
    % solution to the interior points - the boundary points were set up
    % previously
    usol(IntPoints,IntPoints) = reshape(U,N,N);
    
    %Error(i) = norm(usol-Uhat,2)/norm(Uhat,2);
    Error(i) = norm(usol-Uhat,Inf)/norm(Uhat,Inf);
end
% % Finding the order of accuracy using least squares
% 
XX = [ ones(length(h),1) , log(h') ];
bb = log(Error');
%v = [K p]
v = (XX'*XX)\XX'*bb;
disp(v)
disp(Error)


figure(1)
clf
loglog(h,Error,'r*')
xlabel('Step size h')
ylabel('Error in solution given step size h')
title('Error vs. Step size')

figure(2);
mesh(x,y,usol);
title('Numerical solution');

figure(3);
mesh(x,y,Uhat);
title('True solution');