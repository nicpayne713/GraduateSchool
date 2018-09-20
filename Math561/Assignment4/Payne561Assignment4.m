% Nicholas Payne
% Math 561 Assignment 4
% Due November 13
%%
% Problem 4
% Solve the ODE, y'(t) = y(t)/t - t^2 / y^2(t)
% With the initial condition that y(1) = 1
% on the interval [1 , 1.7]
%
% a. Euler 
%
% y_n+1 = y_n + h f(x,y)
%
t = 1; %initial t value
t_final = 1.7; % end of interval
h = 0.05; % stepsize
y = 1; % initial condition
while t <= t_final;
    y= y + h * ( y/t - (t/y)^2);
    t = t+h;
end
y


%%
% b. classical 4th order Runge-Kutta
%
clear all
clc
t = 1; %initial t value
t_final = 1.7; % end of interval

y = 1; % initial condition
h = 0.05;
f(t , y) = @(t,y) y./t - (t.^2./y.^2);
while t < t_final;
k_1 = f(t          , y              );
k_2 = f(t + (1/2)*h, y + (1/2)*h*k_1);
k_3 = f(t + (1/2)*h, y + (1/2)*h*k_2);
k_4 = f(t +       h, y +       h*k_3);
phi = (k_1 + 2*k_2 + 2*k_3 + k_4) / 6;


    y = y + h * phi;
    t = t+h;
end
y
%%

clc;                                               
clear all;

h=0.05;                      % step size
t = 1:h:1.7;                 % Calculates upto y(3)
y = zeros(1,length(t)); 
y(1) = 1;                     % initial condition
f_ty = @(t,y) y/t - t^2/y^2;  

for i=1:(length(t)-1)         % calculation loop
    k_1 = f_ty(t(i),y(i));
    k_2 = f_ty(t(i)+0.5*h,y(i)+0.5*h*k_1);
    k_3 = f_ty((t(i)+0.5*h),(y(i)+0.5*h*k_2));
    k_4 = f_ty((t(i)+h),(y(i)+k_3*h));
    phi = (1/6)*(k_1+2*k_2+2*k_3+k_4);

    y(i+1) = y(i) + phi*h;  
end
y(i+1)
% c. ode45

%%
% Problem 5

