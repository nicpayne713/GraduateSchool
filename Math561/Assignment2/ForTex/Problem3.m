%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%
fs = @(x) log(x); %the function that we want to interpolate

p = @(x) log(10) *((x - 11)* (x - 12)) / 2 ...
+ log(11)* ((x - 10) *(x - 12)) / -1 + log(12)*(x - 10) *(x - 11)/2;
% p is constructed using Lagrange polynomial interpolation method

figure 
fplot(fs, [10,12], 'g')
hold on
fplot(p, [10,12], '--r')
title('Problem 3 - interpolation vs. ln(x)')
legend('ln(x)', 'interpolation')
hold off

figure 
fplot(fs, [11.0999 11.1001], 'g')
hold on
fplot(p, [11.0999 11.1001], 'r')
legend('ln(x)', 'interpolation')
title('Problem 3 - zoomed in at x=11.1')
% error: e_n (x) = \omega(x)/(n+1)! * f^(n+1)(\xi)
% To find the minimum and maximum bounds, we take the 3rd derivative of the
% function (since we used 3 points for interpolation) and evaluate the error
% evalute the error function at 10 and 12 for the lower and upper bounds.
syms x
f = log(x); % approximated function
f3 = diff(diff(diff(f))); % that function's 3rd derivative
poly2sym(f3, 'x');
f3 = matlabFunction(f3); % make f3 a function in matlab
omega = @(x) (x-10)*(x-11)*(x-12); % the omega function from the error formula
r = omega(11.1)/6; % evaluate omega at 11.1. 6 is (n+1)!
m = f3(10)*r; % f3*r finishes the error formula, this is the lower bound
M = f3(12)*r; % upper bound
% after running code, here are the outputs for the upper and lower bounds
%
% m =  -3.299999999999989e-05 
% M =  -1.909722222222215e-05
%
