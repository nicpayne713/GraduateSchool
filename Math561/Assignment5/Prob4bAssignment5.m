%% Problem 4
%
% PART B
%
clc
clear all
close all
syms k1 k2 h yn y lambda;
% Apply the method to the test problem y' = h*lambda
[k1 , k2] = solve(k1 == lambda*(yn + (1/4)*h*k1 - (1/4)*h*k2), ...
    k2 == lambda*(yn + (1/4)*h*k1 + (5/12)*h*k2), k1, k2)

[y] = solve(y == yn + h*((1/4)*k1 + (3/4)*k2), y)
%to find the reduced characteristic equation, create a function
%to represent y_{n+1} = y_n+h*phi and set h and lambda = 0
yn1 = matlabFunction(y) 
% when Matlab creates a function from a symbolic expression it lists the
% inputs in alphabetical order, so yn1 is a function of h,l,yn
a=yn1(0,0,yn) % yields that y_{n+1}=y_n which tells us that k1 and k2 were 
% simply multiples of yn, so then the reduced characteristic equation is 
% r -1 = 0. Therefore the root is 1 so the method is stable.