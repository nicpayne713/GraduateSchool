%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Find the polynomial of degree 4 that solves the generalized interpolation
% problem: p(0) = 2, p(1) = -4,p(2) = 44, p'(0) = -9, and p'(1) = 4
% PART a) Set up and solve a system of equations
% p(x) = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
% p'(x) = 4*a4*x^3 + 3*a3*x^2 + 2*a2*x + a1

% Set up a matrix, A, to sovle the system
A = [1 0 0 0 0
    1 1 1 1 1
    1 2 4 8 16
    0 1 0 0 0
    0 1 2 3 4]; 
B = [2; -4; 44; -9; 4];

x = A\B;

% These are the coefficient values given in the output of this program
% a0 =  2
% a1 = -9
% a2 =  1
% a3 = -3
% a4 =  5



% PART b) Use a divided difference table
% xi yi/y'i
% 0    2   -9    3    7    5
% 0    2   -6    10   17
% 1   -4    4    44   
% 1   -4    48
% 2    44
%
% So we have that the coefficients of the Newton Polynomial are b_i where
% b_0 = 2, b_1 = -9, b_2 = 3, b_3 = 7, b_4 = 5
%
% Then the Newton Polynomial is 2 - 9x + 3x^2 + 7x^2(x-1) + 5x^2(x-1)^2
% Using a CAS we see this simplifies to give us the 4th degree polynomial
% 5x^4 - 3x^3 + x^2 - 9x +2
% Notice that these coefficients match the results from the system of
% equations given above.
