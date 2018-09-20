%517 Assignment 1 Part 3
% Finding the 4th order accurate CFD for the second derivative
syms h du d2u d3u d4u 
syms du5mu ux0 ux1 ux2 ux3 ux4
O = [1            1        1  1         1;
    -2*h         -h        0  h         2*h;
    2*h^2        .5*h^2    0 .5*h^2     2*h^2;
    -(8/6)*h^3  -(1/6)*h^3 0 (1/6)*h^3  (8/6)*h^3;
    (16/24)*h^4 (1/24)*h^4 0 (1/24)*h^4 (16/24)*h^4];
v = [0 0 1 0 0]';
omega = O\v;
disp(omega)
% apply the IVT so collapse the error term to only depend on 1 c(x)
tau = h^5/120 * du5mu/(-32*omega(1)-omega(2)+omega(4)+32*omega(5));
% get a simplified expression for the error based on the omegas
disp(tau)
% get a simplified expression for u''
d2u = simplify(omega(1)*ux3 + omega(2)*ux4 + omega(3)*ux0 + omega(4)*ux1 + omega(5)*ux2);
disp(d2u)

% Construct the method for this BVP based on d2u from above
h = input('Enter 4 step sizes: ');
    N = 1./h - 1 ;
    RelError = zeros(length(N),max(N)+2);
    Exact = zeros(length(N),max(N)+2);
    Approx = zeros(length(N),max(N)+2);
    AbsError = zeros(length(N),max(N)+2);
    % This for loop creates a new matrix for each input of h vector, uses
    % that to come up with the approximation then stores the errors and
    % exact values in the matrices denoted above 'RelError' , 'Exact', and
    % 'Approx'. The plots do not use the matrices exactly since the first
    % few rows of those matrices are padded with tons of zeros. That's
    % addressed after the loop, they are just for storing values.
    for i = 1:4
    e = ones(N(i)+1,1);
    A = spdiags([1*e, -16*e, (12*h(i)^2+30)*e, -16*e, 1*e], -2:2, N(i)+1 , N(i)+1);
    % Need to modify the matrix due to the periodic conditions
    A(1,:) = [12*h(i)^2+30, -16, 1, zeros(1,N(i)-4), 1, -16];
    A(2,:) = [-16, 12*h(i)^2+30, -16, 1, zeros(1,N(i)-4), 1];
    A(N(i),:) = [1, zeros(1,N(i)-4), 1, -16, 12*h(i)^2+30, -16];
    A(N(i)+1,:) = [-16, 1, zeros(1,N(i)-4), 1, -16, 12*h(i)^2+30];
    disp(min(eig(full(A))))
    % The exact solution is U(x) = sin(4pix)/(1+16pi^2)
    x = linspace(0,1-h(i),N(i)+1); %don't include 1, f(1) will be appended later
    t = linspace(0,1,N(i)+2);
    xx = 4*pi.*x;
    tt = 4*pi.*t;
    d = 1+16*pi^2;
    Uhat = sin(tt)./d; % exact solution gets evaluated at 1, but the approximations
    % don't because of the periodic conditions, so those vectors are
    % just appended a little further down
    f = sin(xx);
    f = 12*h(i)^2.*f;
    u = A\f';
    disp(norm(A*(u'-sin(xx)./d)',2)) % this is - tau
    u = [u;u(1)]; f = [f,f(1)];
    zeropad = zeros(1,max(N)+2 - length(u'));
    % Now the loop stores the calculated in the respective rows of the
    % respective matrices for plotting purposes. The zero pad allows the
    % calculations done with n(i) < max(n) to be put into the matrices
    % because the dimensions have to line up exactly.
    RelError(i,:) = [((u' - Uhat)./u') , zeropad];
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
    xx1 = linspace(0,1,N(1)+2);
    xx2 = linspace(0,1,N(2)+2);
    xx3 = linspace(0,1,N(3)+2);
    xx4 = linspace(0,1,N(4)+2);
    
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
    
%%%%%%%%% Verify 4th order accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
a = (X'*X)\X'*b;
disp(a)
    
    