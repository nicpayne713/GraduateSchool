% Solve the ODE
% y' = y^2(1-y)
% y(0) = 10^-4
% on the interval [0,20000] using both ode45 and ode23s with a requested
% accuracy of 10^-4. Count the number of steps taken by each subroutine on
% the interval [0,10000],and [10000,20000].

dy = @(t,y) y.^2 .* (1-y);
y0 = .0001;
tspan = [0 , 20000];
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4)
[T , Y]=ode45(dy , tspan , y0, options );

counter = sum(T<=10000)-1 % subtract 1 for the initial condition

counter2 = length(T) - counter -1 % steps between 10000 and 20000

x = linspace(10000 , 20000 , length(T));
yy = spline(T , Y , x);
figure
plot(T , yy , '*', T, yy, '-')
title('6 ODE45 solution')
xlabel('Points at which ODE45 reevaluates the ODE')
ylabel('Outputs based on ODE45 step sizes')

[T1 , Y1]=ode23s(dy , tspan , y0, options );

counter23s = sum(T1<=10000)-1 % subtract 1 for the initial condition

counter223s = sum(T1>= 10000) % steps between 10000 and 20000
x1 = linspace(10000 , 20000 , length(T1));
yy1 = spline(T1 , Y1 , x1);
figure
plot(T1 , yy1 , '*', T1, yy1, '-')
title('6 ODE23s solution')
xlabel('Points at which ODE23s reevaluates the ODE')
ylabel('Outputs based on ODE23s step sizes')
