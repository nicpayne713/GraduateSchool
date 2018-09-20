% Problem 5.6.2
clear all
clc
% The coefficient matrix for our system 
A = [105, 140, 210, 420
    84, 105, 140, 210
    70, 85, 105, 140
    60, 70, 84, 105];
%
% Target vector
B = -[840
    1260
    2100
    3780];
%
% Want to solve Ap = B
% The solution vector, our polynomial p(t)
q = A\B;
%
p = poly2sym((q'),'t')
% Output was
% p =
% - (23800*t^3)/23 + (18228*t^2)/23 + (2640*t)/23 - 1492/23
