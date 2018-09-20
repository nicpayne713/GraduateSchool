% c. ode45
%
odefun = @(t, y)  y./t - t.^2./y.^2;  
tspan = [1 , 1.7];
%x = linspace(1 , 1.7);
y0 = 1;
[T , Y] = ode45(odefun , tspan , y0);
%x = linspace(1 , 1.7 , length(T));
x = 1:0.05:1.7;
yy = spline(T , Y , x);
plot(x , yy , '*', x, yy, '-')
title('4c ODE45 solution')
xlabel('Points at which ODE45 reevaluates the ODE')
ylabel('Outputs based on ODE45 step sizes')
steps = length(T) - 1 % need to subtract 1 for the initial value

