%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% f(x) = 1   : int(0,1) x^(1/3) = 3/4 = w0+w1+w2
% f(x) = x   : int(0,1) x^(4/3) = 3/7 = .5w1+w2
% f(x) = x^2 : int(0,1) x^(7/3) = 3/10 = .25w1+w2
A = [1,1,1 ; 0,.5,1; 0,.25,1]; % weight coefficient matrix
v = [3/4 ; 3/7 ; 3/10]; % exact integrals for f(x)'s
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