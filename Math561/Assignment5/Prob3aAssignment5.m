% Assignment 5
% Problem 3
% PART a
%
% solve the given ode with a 4th order Adams-Bashforth-Moulten method using
% Runge-Kutta startup
% stepsize h = 0.05
% interval [a,b] = [1 , 1.7]
%
% first use 4th order RK to get some initial points
clc
%close all
clear all
h=0.05;                      % step size
h24 = h/24;
t = 1:h:1.7;                
y = zeros(1,length(t));     % vector for the corrector

y(1) = 1;                     % initial condition

f = @(t,y) y/t - t^2/y^2;  

for i=1:3         % calculation loop
    k1 = f(t(i),y(i));
    k2 = f(t(i)+0.5*h,y(i)+0.5*h*k1);
    k3 = f((t(i)+0.5*h),(y(i)+0.5*h*k2));
    k4 = f((t(i)+h),(y(i)+k3*h));
    phi = (1/6)*(k1+2*k2+2*k3+k4);
     
    y(i+1) = y(i) + phi*h;
end
% vector to store function evaluations
fval = zeros(1,length(t));

fval(1) = f(t(1),y(1));
fval(2) = f(t(2),y(2));
fval(3) = f(t(3),y(3));
fval(4) = f(t(4),y(4));

for i = 4:length(t)-1;
    %Predictor
    y(i+1)= y(i) + ...
        (55 .* fval(i) - 59 .* fval(i-1) + ...
        37 .* fval(i-2) - 9 .* fval(i-3))...
        .* h24;
    %Evaluate
    fval(i+1) = f(t(i+1), y(i+1));
    %Corrector
    y(i+1) = y(i) + ...
        (9 .* fval(i+1) + 19 .* fval(i) - ...
        5 .* fval(i-1) + fval(i-2))...
        .* h24;
     %Evaluate
    fval(i+1) = f(t(i+1), y(i+1));
end
% now calculate the exact values
true = inline('t .* sign(1 - 3*log(t)) .* abs(1 - 3 * log(t)).^(1/3)');

figure
plot(t, true(t), '*-r',  t, (y), 'go-')
title('3a Solution using exact solution and solution form 4th order ABM method')
xlabel('Points on interval 1 to 1.7 with h= 0.05')
ylabel('Outputs of exact solution vs. ABM method')
legend('true value' , 'approximation')
%
figure
plot(t , true(t)-(y), '--r')
title('3a Difference between exact solution and 4th order ABM method')
xlabel('Points on interval 1 to 1.7 with h - 0.05')
ylabel('Error between method and exact value')
% print out the results
disp(' '); disp('Solution Values'); disp(' '); disp(' t     ABM              exact');...
    disp(' ================================================================');...
    for i = 1:length(t) disp(sprintf(' %4.2f%15.10f%15.10f%15.10f%15.10f',... 
        t(i),y(i),true(t(i)))); end
disp(' '); disp('Errors'); disp(' '); disp(' t       ABM ');...
    disp(' =================================================');...
    for i = 1:length(t) disp(sprintf(' %4.2f%15.10f%15.10f%15.10f',t(i),y(i)-true(t(i)) ... 
       )); end
