%Problem 4
% Solve the ODE, y'(t) = y(t)/t - t^2 / y^2(t)
% With the initial condition that y(1) = 1
% on the interval [1 , 1.7]
%
% d. by hand
%
% y_n+1 = y_n + h f(x,y)
h=0.05;                      % step size
t = 1:h:1.7;               % points to evaluate at 

true = @(t) t.*((1-3.*log(t)).^(1/3));
[Y] = true(t)

figure
plot(t , Y , '*', t, Y, '-r')
title('4d Solution y(t) using exact solution')
xlabel('Points used in parts a and b')
ylabel('Outputs of exact solution')