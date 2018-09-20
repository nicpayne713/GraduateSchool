%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% create a coefficient matrix, W, for the weights and solve Wv=a
% where v is the vector [w0 w1 w2 w3]' and a is the target vector 
% a=[0 1/h 0 0]'
% We want our equation for f'(x) to only equal f'(x) so we kill off the
% other terms by setting them equal to 0 and use the 1/h factor to get rid
% of the hf' in the equation.
syms h
W = [1 1 1 1
0 1 2 4
0 .5 2 8
0 1/6 4/3 32/3];

a = [0
8 % Since our step size is h = 1/8, to solve for the weights directly I 
% set h=1/8 so that the vector v will contain the numerical weight values
% If we were taking different initial values, then we would simply change
% this value to be 1/h again where h is still the step size.
0
0];

v = W\a;
w0 = v(1);
w1 = v(2);
w2 = v(3);
w3 = v(4);

f = @(x) exp(x);
% f' = f
true = f(0)
approx = w0*f(0) + w1*f(1/8) + w2*f(1/4) + w3*f(1/2)
trueerror = true - approx