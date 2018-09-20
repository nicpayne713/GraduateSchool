% Math 517
% Assignment 1
% Part 2 (#4-6)
% First I create a sparse matrix A to solve the system AU = f for the given
% BVP
clc 
clear 

    h = input('Enter 4 step sizes: ');
    N = 10./h - 1 ;
    RelError = zeros(length(N),max(N)+2);
    Exact = zeros(length(N),max(N)+2);
    Approx = zeros(length(N),max(N)+2);
    AbsError = zeros(length(N),max(N)+2);
    % This for loop creates a new matrix for each input of n vector, uses
    % that to come up with the approximation then stores the errors and
    % exact values in the matrices denoted above 'RelError' , 'Exact', and
    % 'Approx'. The plots do not use the matrices exactly since the first
    % few rows of those matrices are padded with tons of zeros. That's
    % addressed after the loop, they are just for storing values.
for i = 1:4
   % h(i) = 10/(N(i)+1);
    % This is the 2nd order finite difference method for the BVP
    e = ones(N(i)+2,1);
    A = spdiags([e (h(i)^2-2)*e e], -1:1, N(i)+2 , N(i)+2);
    A(1,1) = h(i)^2-2*h(i)-2;
    A(1,2) = 2;
    A(N(i)+2,N(i)+1) = 2;
    A(N(i)+2,N(i)+2) = h(i)^2-2*h(i)-2;
% B 
    %x = 0:h(i):10;
    x = linspace(0,10,N(i)+2 );
% The exact solution is:
    c = exp(10)/(2*cos(10));
    Uhat = c.*sin(x) + c.*cos(x) - .5.*exp(x) ;
    f = (-1).*exp(x);
    f = h(i)^2.*f; % the matrix A doesn't include the h^2 factor - it was moved to the RHS
    u = A\f';
    zeropad = zeros(1,max(N)+2 - length(u'));
    % Now the loop stores the calculated in the respective rows of the
    % respective matrices for plotting purposes. The zero pad allows the
    % calculations done with n(i) < max(n) to be put into the matrices
    % because the dimensions have to line up exactly.
    RelError(i,:) = [((u'-Uhat)./u') , zeropad];
    Exact(i,:) = [Uhat , zeropad];
    Approx(i,:) = [u' , zeropad];
    AbsError(i,:) = abs([((u'-Uhat)) , zeropad]);
end
    x1 = ['mesh points with  h = ', num2str(1/(N(1)+1))];
    y1 = ['error'];
    x2 = ['mesh points with  h = ', num2str(1/(N(2)+1))];
    y2 = ['error'];
    x3 = ['mesh points with  h = ', num2str(1/(N(3)+1))];
    y3 = ['error'];
    x4 = ['mesh points with  h = ', num2str(1/(N(4)+1))];
    y4 = ['error'];
    % These grab the necessary values of the Relative Error, Exact solution
    % , and the approximate solution. If you plot the whole row of each
    % matrix for each value, then the 1st, 2nd, and 3rd rows have a ton of
    % zeros which mess up the results, so this makes it so that it graphs
    % only the calculations and not the zeros used to pad the matrix when
    % it was being populated in the loop
    RelError1 = RelError(1,(1:N(1)+2));
    RelError2 = RelError(2,(1:N(2)+2));
    RelError3 = RelError(3,(1:N(3)+2));
    RelError4 = RelError(4,(1:N(4)+2));
    Exact1 = Exact(1,(1:N(1)+2));
    Exact2 = Exact(2,(1:N(2)+2));
    Exact3 = Exact(3,(1:N(3)+2));
    Exact4 = Exact(4,(1:N(4)+2));
    Approx1 = Approx(1,(1:N(1)+2));
    Approx2 = Approx(2,(1:N(2)+2));
    Approx3 = Approx(3,(1:N(3)+2));
    Approx4 = Approx(4,(1:N(4)+2));
    xx1 = linspace(0,10,N(1)+2);
    xx2 = linspace(0,10,N(2)+2);
    xx3 = linspace(0,10,N(3)+2);
    xx4 = linspace(0,10,N(4)+2);
figure(1)
    clf
    subplot(2,2,1)
    plot(xx1,RelError1,'r-')
    title('Relative Error between exact solution on the mesh and the approximation')
    xlabel(x1)
    ylabel(y1)
    
    subplot(2,2,2)
    plot(xx2,RelError2,'r-')
    title('Relative Error between exact solution on the mesh and the approximation')
    xlabel(x2)
    ylabel(y2)
    
    subplot(2,2,3)
    plot(xx3,RelError3,'r-')
    title('Relative Error between exact solution on the mesh and the approximation')
    xlabel(x3)
    ylabel(y3) 
    
    subplot(2,2,4)
    plot(xx4,RelError4,'r-')
    title('Relative Error between exact solution on the mesh and the approximation')
    xlabel(x4)
    ylabel(y4)
    
figure(2)
    clf
    subplot(2,2,1)
    plot(xx1,Exact1,'r-',xx1,Approx1, 'b-')
    title('Exact vs. Approximation given spacing noted on x axis')
    xlabel(x1)
    ylabel('Outputs of solutions')
    
    subplot(2,2,2)
    plot(xx2,Exact2,'r-',xx2,Approx2, 'b-')
    title('Exact vs. Approximation given spacing noted on x axis')
    xlabel(x2)
    ylabel('Outputs of solutions')
    
    subplot(2,2,3)
    plot(xx3,Exact3,'r-',xx3,Approx3,'b-')
    title('Exact vs. Approximation given spacing noted on x axis')
    xlabel(x3)
    ylabel('Outputs of solutions')
    
    subplot(2,2,4)
    plot(xx4,Exact4,'r-',xx4,Approx4,'b-')
    title('Exact vs. Approximation given spacing noted on x axis')
    xlabel(x4)
    ylabel('Outputs of solutions')
    
%%%%%%%%% Verify 2nd order accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%
AbsError1 = AbsError(1,(1:N(1)+2));
AbsError2 = AbsError(2,(1:N(2)+2));
AbsError3 = AbsError(3,(1:N(3)+2));
AbsError4 = AbsError(4,(1:N(4)+2));
Error = [norm(AbsError1,Inf),norm(AbsError2,Inf),norm(AbsError3,Inf),norm(AbsError4,Inf)];


figure(3)
scatter(h,Error)
lsline

X = [ ones(4,1) , log(h') ];
b = log(Error');
%v = [K p]
v = (X'*X)\X'*b;
disp(v)