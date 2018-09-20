function [e_z] = raytrace(theta)
% This function takes an angle in degrees as its input
% Snell's Law
global A c cPrime dz2dx

%A = cos(degtorad(theta0))/z0;
x=linspace(0,145824);
% ODE
%dz2dx = @(t,y) [y(2);-cPrime(y(1))./((A(2000,5.4)).^2.*c(y(1)).^3)];
e_z=[]; % initialize matrix to hold all values in

if nargin==0 %checks for argument
    warning('"raytrace(x)" must contain an argument. Else, answer given on linspace.')
    return
end
format long
% Initial conditions
y = zeros(1,2);
y(1) = 2000;% this is the initial condition from part A, I'm afraid that 
            %the initial condition actually depends on theta
y(2) = tan(degtorad(theta));

    % receiver is located at x = 24 miles = 145824 ft.
    % receiver is located at depth of z = 3000 ft.
    [~ , Y] = ode45(dz2dx, x, y);
    e_z = Y(end,1)-3000; % find the distance that the ray missed the receiver
    
end
