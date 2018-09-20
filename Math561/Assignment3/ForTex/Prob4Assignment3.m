%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% find the roots of p_2(x) to use in this script. Using a CAS, the exact
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
% now we apply the approximation to the requested function
f = @(x) exp(-x);
true = gammainc(1,4/3) * gamma(4/3)
approx = w0*f(X(2,1)) + w1*f(X(2,2))