% Math 517 Assignment 2 Part 1
clear all
close all
% Construct sparse coefficient matrix A and appropriate right-hand side
% vector F for the descretized version utilizing the 5-point Laplacian on a
% uniform mesh of Poisson's equation in 2D 
% First I will set up the grid
a = 0;
b = 1;
% The grid will be based on the desired step size for the problem
% For the purposes of doing the convergence study, h is a vector and the
% script loops through the problem for how ever many h's there are to try.
h = [1/10,1/100, 1/1000];
Error = zeros(1,length(h));
for i = 1:length(h)
    N = (1/h(i))-1;
    % On page 68 of LeVeque, he includes a MATLAB code which will help
    % construct the sparse matrix A for Laplacian(u) = f using the
    % kronecker tensor product.
    I = speye(N);
    e = ones(N,1);
    T = spdiags([e, -4*e, e],-1:1, N,N);
    S = spdiags([e e], [-1 1], N,N);
    A = (kron(I,T) + kron(S,I));
    A = -A; % The sign for A gets switched because this problem has the
            % minus Laplacian
    
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
    
    % f needs to be evaluated at all the interior points, so that's done here
    RHS = h(i)^2 * f(Xinterior,Yinterior);
    
    % The boundary conditions are taken care of below
    % Uhat is the exact solution sampled on the whole mesh
    Uhat = exp(X+0.5.*Y);
    
    % Now we set the boundary conditions, since we have a known solution around
    % the boundary we just need a full array and we'll use the boundary
    % values only
    usol = Uhat;
    
    % Now we modify the RHS so that the boundary terms are included
    RHS(:,1) = RHS(:,1) - usol(IntPoints,1);
    RHS(:,N) = RHS(:,N) - usol(IntPoints,N+2);
    RHS(1,:) = RHS(1,:) - usol(1,IntPoints);
    RHS(N,:) = RHS(N,:) - usol(N+2,IntPoints);
    
    % Now I can take the grid created with the meshgrid function and stretch it
    % out into an array for the RHS of the system
    F = -reshape(RHS,N*N,1);
    
    % Now to solve the system using the backslash
    U = A\F;
    
    % The U vector gets put into the usol grid as the
    % solution to the interior points - the boundary points were set up
    % previously
    usol(IntPoints,IntPoints) = reshape(U,N,N);
   
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