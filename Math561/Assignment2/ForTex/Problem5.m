%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
% PART A
%
x5 = [-2 0 2 3 4 5]; % inputs
y5 = [4 0 -4 -30 -40 -50]'; % outputs of function
V = vander(x5); % declare Vandermonde matrix to find polynomial
c = V\y5; % vector with coefficients of interpolating polynomial
p = poly2sym(c, 'x'); 
a = -2:.1:5; % domain of the interpolation
polynomial = matlabFunction(p); % fa is the polynomial interpolation
%
% PART B
%
s = spline(x5,y5,a);
%
% PART C
%
% for clamped spline, store end conditions into new vector 'e'
e = [-2 y5' 10]';
clamped = spline(x5,e, a);
%
% PART D
%
chip = pchip(x5, y5, a);

figure 
hold on
scatter(x5,y5,'o', 'c')
plot(a, polynomial(a), ':r')
plot(a, s, '--g')
plot(a, clamped, ':')
plot(a, chip, '--')
legend('data', 'polynomial', 'spline', 'clamped spline', 'pchip')
title('Problem 5')
%