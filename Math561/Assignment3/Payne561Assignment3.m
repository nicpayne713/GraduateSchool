% Nicholas Payne
% Math 561 Assignment 3
% Due October 23, 2014
clear all
close all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% After using Taylor expansion on points 0,h,2h,4h and collecting the like
% terms, we are left with the formula: f'(x) = (w0 + w1 + w2 + w3)f 
% +(w1 + 2w2 + 4w3)hf'
% +(1/2w1 + 2w2 + 8w3)h^2f''
% +(1/6w1 + 4/3w2 + 32/3w3)h^3f'''
% +(1/24w1 + 2/3w2 + 32/3w3)h^4f'''' <-- This is the error term
% So we create a coefficient matrix, W, for the weights and solve Wv=h
% where v is the vector [w0 w1 w2 w3]' and h is the target vector 
% a=[0 1/h 0 0]'
% We want our equation for f'(x) to only equal f'(x) so we kill off the
% other terms by setting them equal to 0 and use the 1/h factor to get rid
% of the hf' in the equation.


 syms h
W = [1 1 1 1
0 1 2 4
0 .5 2 8
0 1/6 4/3 32/3];

a = [0
8 % Since our step size is h = 1/8, to solve for the weights directly I 
% set h=1/8 so that the vector v will contain the numerical weight values
% If we were taking difference initial values, then we would simply change
% this value to be 1/h again where h is still the step size.
0
0];

v = W\a;
w0 = v(1);
w1 = v(2);
w2 = v(3);
w3 = v(4);

f = @(x) exp(x);
% f' = f
true = f(0)

approx = w0*f(0) + w1*f(1/8) + w2*f(1/4) + w3*f(1/2)

% f^(4) = f
error = (w0/24 +2*w2/3 + w3*32/3)*(1/8)^4

% uppbound on error


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PART a)
f = @(x) exp(x)/x;
truea = expint(-1)-expint(-2);
ra = romberg(f, 1, 2, 5);
Ea = ra-truea;
for i = 1:5;
    for j = 1:5;
        if i<j;
            Ea(i,j)=0;
        end
    end
end
ra
Ea
%
% PART b)
g = @(x) sqrt(1-x.^2);
trueb = integral(g, 0,1);
rb = romberg(g, 0,1,5);
Eb = rb-trueb;
for i = 1:5;
    for j = 1:5;
        if i<j;
            Eb(i,j)=0;
        end
    end
end
rb
Eb
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% f(x) = 1   : int(0,1) x^(1/3) = 3/4 = w0+w1+w2
% f(x) = x   : int(0,1) x^(4/3) = 3/7 = .5w1+w2
% f(x) = x^2 : int(0,1) x^(7/3) = 3/10 = .25w1+w2
A = [1,1,1 ; 0,.5,1; 0,.25,1]; % weight coefficient matrix
v = [.75 ; 3/7 ; 3/10]; % exact integrals for f(x)'s
w = A\v; % finds w0,w1,w2
weights = rats(w) % returns the weights as fractions
%
% if f(x) = e^(-x)
% f(0) = 1 , f(.5) = e^(-.5) , f(1) = 1/e
%
% int(0,1) x^(1/3)e^(-x) ~ w0*f(0) + w1*f(.5) + w2*f(1)
f = @(x) exp(-x);
Approx = w(1)*f(0) + w(2)*f(.5) + w(3)*f(1)
true = gammainc(1,4/3) * gamma(4/3)
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% P is the coefficient polynomial for finding p_2(x)
P = [3/4 3/7
   3/7 3/10 ];
% b is the vector containing the constant terms of the inner products of
% the weight function, p_2(x), and f(x)
b = [-3/10;-3/13];
% c is the vector with coefficients a and b for p_2(x) = a + bx + x^2
c = P\b;
% reformat c so that it is the correct polynomial, p_2(x)
c = [1,c(2),c(1)];
% find the foots of p_2(x) to use in this script. Using a CAS, the exact
% roots were found and used for the matrix X, used to solve the system of
% equations involving the weight functions.
r = roots(c);
X = [1 1
    (35-3*sqrt(35))/65 (35+3*sqrt(35))/65];
a = [3/4;3/7];
% creates a vector with w0 and w1
weights = X\a;
w0 = weights(1);
w1 = weights(2);
%
% now we apply the approximation to the requested function
f = @(x) exp(-x);
true = gammainc(1,4/3) * gamma(4/3)
approx = w0*f(X(2,1)) + w1*f(X(2,2))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PART a)
%
f5 = @(x) exp(2.*x)./((1+x.^2).^2);
q = quad(f5, 0,3)
ql = quadl(f5,0,3)
qgk = quadgk(f5,0,3)
int = integral(f5, 0,3)
%
% PART b)
%


%
% PART c)
%
fc = @(x,y) 2.*y.*sin(x)+(cos(x)).^2;
c = @(x) sin(x);
d = @(x) cos(x);
quad2d(fc,0,pi/4,c,d)
integral2(fc,0,pi/4,c,d)