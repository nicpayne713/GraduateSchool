%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 7 EXTRA CREDIT %%%%%%%%%%%%%%%%%%%%%%
%
% In MuPad, we use the Pade command to find the Pade approximation of f(x)
% = tan(x) around x=0 witn m=n=2.
% PART A
% Pade(tan(x),0,[2,2]) = 
% x(x^2-15) / 3(2x^2-5)
% The Pade approximation has a singularity at x = +- sqrt(5/2)
t = [(-sqrt(5/2)) :.1: (sqrt(5/2))];
y = [-2:.1:2];
P = @(x) (x.*(x.^2-15)) ./ (3.*(2.*x.^2-5));
f = @(x) tan(x);
P1 = (y.*(y.^2-15)) ./ (3.*(2.*y.^2-5));
figure

plot(y, P(y), 'or')
hold on
plot(y, tan(y),'g')
ylim([-10 10])
legend('Pade', 'tangent')
title('Pade approximation at 0 with n=m=2 of tan(x)')
%
% PART B
% In MuPad, we use the Pade command to find the Pade approximation of f(x)
% = tan(x) around x=0 witn m=n=3.
% 
% The Pade approximation with m=n=3 is exactly the same as m=n=2 