% 561 Final Assignment 
% Problem 4
% PART a
clear all
close all
clc
%global y;
s = linspace(-40,-5);
tspan=[0,1];
y = g(s); %calls function g which uses values of s as the initial guess for y'(0) 
% and uses the given initial condition for y(0) = 4, then solves the system
% of ODE's broken down from y'' = 3/2y^2. g returns the values of 
%y(s) - 1 to calculate by how much the IVP with 's' as the initial guess
%misses the true solution.

% PART b
% plots g(s) for s in [-40,-5] and the 2 solutions
figure
plot(s,y, 'og',s ,y , '-r')
title('error between given initial value for y(1) in BVP , and solution to IVP based on initial guess for yprime(0)=s')
xlabel('initial guess for yprime(0) = s')
ylabel('error between given solution and calculated solution')
% find the indices for which values of s gave close answers to the solution
% to use as the initial guess for fzero in part c
indices = find([0 diff(sign(y))]~=0)
% indices are 13 and 92
x0 = s(12);
x00 = s(91);
% PART c
x1 = fzero(@g,x0);
x2 = fzero(@g,x00);
% x1 and x2 are the best solutions to y(1) ~ 1 and y(1) ~ 1 based on
% initial guess y'(0) = x1, y'(0) = x2.
% PART d
% ODE function
dy2dt = @(t,z) [z(2) ; (3/2).* z(1).^2];
[T,Y1] = ode45(dy2dt,tspan,[4,x1]);
[T,Y2] = ode45(dy2dt,tspan,[4,x2]);
figure
plot(T,Y1(:,1),'co-',T,Y1(:,2),'bo-')
legend('solution to y based on yprime(0)=x1',...
    'solution to yprime based on yprime(0)=x1')
figure
plot(T,Y2(:,1),'*r--',T,Y2(:,2),'*g--')
legend('solution to y based on yprime(0)=x2',...
    'solution to yprime based on yprime(0)=x2')