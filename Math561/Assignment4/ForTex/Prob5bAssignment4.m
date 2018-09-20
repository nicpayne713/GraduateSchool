%%
% Problem 
% Solve the ODE, y''(t) = -y'(t)/t- (t^2-1)/t * y(t) 
% With the initial condition that y(1) = 1, y'(1) = 0
% on the interval [1 , 2]
% h = 0.1
%
% b. RK4 
clear all
close all
h = 0.1; % stepsize
t = 1:h:2; % interval
y = zeros(1, length(t));
p = zeros(1, length(t));
f1 = @(t , y , p) p;
f2 = @(t , y , p) -p./t - (t.^2-1)./(t.^2)*y;
y(1) = 1;
p(1) = 0;

for i = 1:(length(t)-1);
    k1 = h * f1(t(i)        , y(i)         , p(i)        );
    l1 = h * f2(t(i)        , y(i)         , p(i)        );
    k2 = h * f1(t(i) + .5*h , y(i) + .5*k1 , p(i) + .5*l1);
    l2 = h * f2(t(i) + .5*h , y(i) + .5*k1 , p(i) + .5*l1);
    k3 = h * f1(t(i) + .5*h , y(i) + .5*k2 , p(i) + .5*l2);
    l3 = h * f2(t(i) + .5*h , y(i) + .5*k2 , p(i) + .5*l2);
    k4 = h * f1(t(i) +    h , y(i) +    k3 , p(i) +    l3);
    l4 = h * f2(t(i) +    h , y(i) +    k3 , p(i) +    l3);
    
    y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    p(i+1) = p(i) + (l1 + 2*l2 + 2*l3 + l4)/6;
end
y
p
x = linspace(1,2);
figure  
plot(t , y, '-c', t , y, '*')
title('The Plot of y(t) using the RK4 Method')
xlabel('Interval used for finding solution')
ylabel('Outputs from RK4 method at each step')
figure
plot(t , p , '-r', t , p, '*')
title('The Plot of y`(t) using the RK4 Method')
xlabel('Interval used for finding solution')
ylabel('Outputs from RK4 method at each step')