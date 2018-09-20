%%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f(x)=xln(x) on [1,3] ----> f(x)=(x+2)(ln(x+2)) on [-1,1]
%Chebyschev Polynomials: T0 = 1, T1 = x, T2 = 2x^2 - 1, T3 = 4x^3 -3x, and
%T4 = 8x^4 -8x^2 +1
%Use the Chebyshev Polynomials to find the best least-squares
%approximation by a polynomial of degree 4  of f. Print those coefficients.
%Note that Ti are not normalized.
%Note that approximating xln(x) on [1 3] is equivilant to approximating
%(x+2)(ln(x+2)) on [-1 1].

f = @(x) (x+2).*(log(x+2)); %The function we are approximating 
w = @(x) 1/((1-x^2)^.5); %Weight function for Cheb. Poly's of 1st kind

%Stores Chebyshev polynomial coefficients as row vectors in matrix T via a
%script which prints out the coefficients for Cheb. Poly's of 1st kind.
T = t_polynomial_coefficients(4); 

%Stores each Chebyshev polynomial into its own vector leading with the 4th
%degree term
T0 = poly2sym(wrev(T(1,:)),'x');
T1 = poly2sym(wrev(T(2,:)),'x');
T2 = poly2sym(wrev(T(3,:)),'x');
T3 = poly2sym(wrev(T(4,:)),'x');
T4 = poly2sym(wrev(T(5,:)),'x');

%To find the constants for each Chebyshev polynomial we call the
%integration operator from MuPad to compute the inner-products and then
%divide by the normalization factor for each Ti
c0 = vpa(int(f*T0*w, -1, 1),10)/(pi); 
c1 = vpa(int(f*T1*w, -1, 1),10)/(pi/2);
c2 = vpa(int(f*T2*w, -1, 1),10)/(pi/2);
c3 = vpa(int(f*T3*w, -1, 1),10)/(pi/2);
c4 = vpa(int(f*T4*w, -1, 1),10)/(pi/2);

% Creates a matrix A with each scaled polynomial
A = [ c0*T0
    c1*T1
    c2*T2
    c3*T3
    c4*T4 ];

S = sum(A); %Sums the columns of A which are the coefficients of each power of x

P = poly2sym(S, 'x'); %Convets vector to sym expression of the polynomial

v = matlabFunction(P); % Converts the sym expression for P to a function v

figure 
t = linspace(-1,1);
plot(t+2, v(t), 'r') % Graphs my approximation in red
hold on
plot(t+2, f(t), ':c') % Graphs the original function
xlabel('-1 <= x <= 3')
ylabel('x*ln(x) and my approximation')
title('Problem 1 - Degree 4 least-squares approximation of x*ln(x)')
legend('approximation', 'x*ln(x)')
%
% plots the error
y1 = v(t);
y2 = f(t);
figure 
plot(t+2, y1-y2, 'r')
title('Problem 1 - error of approximation on [1,3]')

% Prints out each chebyshev polynomial scaled by its coefficient
A(1);
A(2);
A(3);
A(4);
A(5);

% Prints out the polynomial which comes from summimg the chebyshev
% polynomials

P;

% After running the code these are the coefficient values, included as a
% comment in this script and also included in a separate output file
% Here is the sum of Chebyshev polynomials
% 4.7613055409573337983797358674565/pi + (5.2141299184740205197030604722386*x)/pi
% +(0.41082067310163193629701794407083*(2*x^2 - 1))/pi
% +(0.036243183818908685650539425182615*(- 4*x^3 + 3*x))/pi
% +(0.004819505329342707589113103949785*(8*x^4 - 8*x^2 + 1))/pi
%
%
%
% Here is the polynomial. My polynomial P approximates xln(x) on -1 to 1 so
% to shift it to 1 to 3 we replace x by x-2, so with the help of a CAS
% we get that P(x-2) is:
% (19278021317370832x^4 - 226710538176784036x^3 + 1289133369228065366x^2
% - 391473965981346621x - 690697394020706867) / (500000000000000000pi)
