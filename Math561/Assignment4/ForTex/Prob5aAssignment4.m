%%
% Problem 
% Solve the ODE, y''(t) = -y'(t)/t- (t^2-1)/t * y(t) 
% With the initial condition that y(1) = 1, y'(1) = 0
% on the interval [1 , 2]
% h = 0.1
%
% a. Euler 
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
    
    y(i+1) = y(i) + k1;
    p(i+1) = p(i) + l1;
end
y
p
figure  
plot(t , y, '-c', t , y, '*')
title('The Plot of y(t) using Eulers Method')
xlabel('Interval used for finding solution')
ylabel('Outputs from Eulers method at each step')
figure
plot(t , p , '-r', t , p, '*')
title('The Plot of y`(t) using Eulers Method')
xlabel('Interval used for finding solution')
ylabel('Outputs from Eulers method at each step')