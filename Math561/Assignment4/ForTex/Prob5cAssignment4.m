% Problem 5
% Solve the ODE, y''(t) = -y'(t)/t- (t^2-1)/t * y(t) 
% With the initial condition that y(1) = 1, y'(1) = 0
% on the interval [1 , 2]
% h = 0.1
%
% c. ODE45 
    
f = @(t, y) [y(2);(1-t^2)*y(1)-t*y(2)/(t^2)];
[T , Y] = ode45(f, [1 2] , [1; 0]); 
x = 1:.1:2;

yy = spline(T, Y(:,1) , x);
figure
plot(x, yy, 'o', x, yy, 'r')
title('The solution y`(t) to the ODE')
xlabel('Evaluating ODE45 solution at points from a and b')
ylabel('Outputs of y`(t)')

yyy = spline(T , Y(:,2), x);
figure
plot(x, yyy , 'o', x, yyy, 'r')
title('The solution y(t) to the ODE')
xlabel('Evaluating ODE45 solution at points from a and b')
ylabel('Outputs of y(t)')

steps = length(T) - 1 % need to subtract 1 for the initial value