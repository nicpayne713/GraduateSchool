format long
clear all
clc
close all
% Approximate f(x) = xln(x)on [1,3] using Chebyshev polynomials
% We want to get rid of the singularity in our inner products, so
% we substitute x = cos(theta) and also represent each
% Chebyshev polynomial T_k as cos(k*theta). 
% Also note, that we will shift our function to be over [-1 1] since that's
% where the Chebyshev polynomials are defined
% c0...c4 are the coefficients of the each Chebyshev polynomial that will
% be summed to give the approximating polynomial per the equation on page
% 90 of Gautschi
F0 = @(x) (cos(x)+2).*log(cos(x)+2).*cos(0.*x); % Function to find c0
F1 = @(x) (cos(x)+2).*log(cos(x)+2).*cos(1.*x); % Function to find c1
F2 = @(x) (cos(x)+2).*log(cos(x)+2).*cos(2.*x); % Function to find c2
F3 = @(x) (cos(x)+2).*log(cos(x)+2).*cos(3.*x); % Function to find c3
F4 = @(x) (cos(x)+2).*log(cos(x)+2).*cos(4.*x); % Function to find c4

% The following are the inner products to find the coefficients of the
% Chebyshev polynomials. Before the trig substitution the inner products
% would be taken between f(x) = x*ln(x), t_k(x) (each Chebyshev polynomial
% of the first kind from 0 to 4), and w(x) = 1/sqrt(1+x^2) which is the
% weight function for Chebyshev polynomials of the first kind, from -1 to 1.
%However, after the trig substitution and the shift on the x-axis from 
%[1 3] to [-1 1], the inner products are taken of Fk (defined above)from 0 to pi
%and multiplied by the normalizing factors (1/pi for the constant Chebyshev 
%polynomial, and 2/pi for the others.
c0 = (1/pi)*quad(F0,0,pi); % The normalizing factor for c0 is 1/pi
c1 = (2/pi)*quad(F1,0,pi); % The normalizing factor for c1...c4 is 2/pi
c2 = (2/pi)*quad(F2,0,pi); 
c3 = (2/pi)*quad(F3,0,pi); 
c4 = (2/pi)*quad(F4,0,pi); 

% We symbolically define the Chebyshev polynomials (found on page 87 of our
% book, or with a simple Google search) in order to print out the sum of
% Chebyshev polynomials which approximates x*ln(x) on [1,3]

syms x
T0 = 1;
T1 = x;
T2 = 2*x^2 -1;
T3 = 4*x^3 -3*x;
T4 = 8*x^4 -8*x^2 +1;

% Let T be the sum of Chebyshev polynomials scaled by each coefficient per
% the formula on page 90 of Gautschi

T = c0*T0 +c1*T1 +c2*T2 +c3*T3 +c4*T4;

% However we need to shift the approximation back onto [-1,3] so we sub into
% the equation for x, x-2

T1 = subs(T,x-2);

% To plot the error we simply plot the difference between T1 and x*ln(x)

l = linspace(1,3);
% I create functions t1 and f1 to calculate the error over [1,3]
t1 = matlabFunction(poly2sym(T1,'x')); 
f1 = @(x) x.*log(x);
yt = t1(l);
yf = f1(l);
plot(l,yt-yf, 'r')
title('error between f(x)=x*ln(x) and Chebyshev Approximation')

figure
fplot(t1,[1 3],'c--')
hold on
fplot(f1,[1 3],'b:')
legend('Chebyshev Approximation', 'f(x) = x*ln(x)')
title('Problem 1 - Chebyshev Approximation of f(x) = x*ln(x)')

% Prints out the coefficients from the inner products taken above
format long
ChebyCoeff = [c0
c1
c2
c3
c4]

PolyApprox = sym2poly(T1)
