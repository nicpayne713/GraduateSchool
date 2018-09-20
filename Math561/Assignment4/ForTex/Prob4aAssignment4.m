%%
% Problem 4
% Solve the ODE, y'(t) = y(t)/t - t^2 / y^2(t)
% With the initial condition that y(1) = 1
% on the interval [1 , 1.7]
%
% a. Euler 
%
% y_n+1 = y_n + h f(x,y)
h=0.05;                      % step size
t = 1:h:1.7;                
y = zeros(1,length(t)); 
y(1) = 1;                     % initial condition
f_ty = @(t,y) y/t - t^2/y^2;  

for i=1:(length(t)-1)         % calculation loop
    phi = f_ty(t(i),y(i));
    y(i+1) = y(i) + phi*h;  
end
figure
plot(t , y , '*', t, y, '-r')
title('4a Solution y(t) using Eulers method')
xlabel('Points where the slope of y is reevaluated')
ylabel('Outputs using Eulers method')