%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
clear all
% The Runge Function is given by R(x) = 1/(1+25x^2)
%
% Vector containing all the points that each function will be evaluated at
e = -1:.02:1;
% PART A
% We use a matrix to find the coefficients of the interpolating polynomial
%
% Create a vector with the x-coordinates of the interpolating points
x = -1:.1:1; 
%
V = vander(x);
%
% Define the runge function to get our interpolating points
R = @(x) (1+25.*x.^2).^(-1);
%
% Declare a column vector which holds the y-coordinates for our interpolating
% points
y = R(x)';
%
c = V\y
%
% After running the code, c (a column vector with all the coefficients of
% the interpolating polynomial) is
% c = 1.0e+06 *
%
%           0.260178630971763
%           -0.000000000000002
%           -1.012094874480108
%            0.000000000000007
%            1.639177410848146
%           -0.000000000000010
%           -1.442944963439323
%            0.000000000000007
%            0.757287092774918
%           -0.000000000000003
%           -0.245249283924509
%            0.000000000000001
%            0.049317542912939
%           -0.000000000000000
%           -0.006119219949421
%            0.000000000000000
%            0.000470846226759
%            0.000000000000000
%           -0.000024143479625
%           -0.000000000000000
%            0.000001000000000
%
% v is the interpolating polynomial of degree 20
v = matlabFunction(poly2sym(c,'x'));
%
% PART B
%
% set xt as a vector containing the Chebyshev points
xt = cos((2.*([20:-1:0])+1).*pi./42);
% create a Vandermonde matrix 
Vt = vander(xt);
yt = R(xt)';
% solves the system and sets ct as the column vector containing the
% coefficients of the interpolating polynomial at the Chebyshev points.
ct = Vt\yt;
% converts ct into an inline function to plot
vt = matlabFunction(poly2sym(ct,'x'));
%
% the coefficients of the interpolating function are:
%
% 1.0e+04 *
%
%          0.646655103655755
%         -0.000000000004668
%         -3.420805498370904
%          0.000000000021239
%          7.775445632006652
%         -0.000000000040194
%         -9.930012492409208
%          0.000000000041009
%          7.823630205987262
%         -0.000000000024462
%         -3.933329690078009
%          0.000000000008663
%          1.263561779281532
%         -0.000000000001769
%         -0.253727314541573
%          0.000000000000192
%          0.030662946960456
%         -0.000000000000009
%         -0.002176233358920
%          0.000000000000000
%          0.000100000000000
%
%
% PART C
% use spline interpolation on the 21 points from part a and evaluate the
% spline at the 101 points in vector 'e'
%
yy = spline(x, y, e); 

figure 
% plots each interpolation and its evaluation at the 101 points on [-1 1]
plot(e,v(e),'r', e, vt(e),'g', e,yy,'c')
ylim([-.2 1.2])
legend('interpolation', 'interpolation at cheby points', 'spline')
title('Problem 4')