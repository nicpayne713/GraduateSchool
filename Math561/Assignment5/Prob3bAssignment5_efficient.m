% Assignment 5
% Problem 3
% PART b
% first use 4th order RK to get some initial points
clc
close all
clear all
h = 0.1; % stepsize
h24 = h/24;
t = 1:h:2; % interval
y = zeros(2, length(t));
y(:,1) = [1;0];
f = @(t,z) [z(2,1);
            -z(2,1)./t - (t.^2-1)./(t.^2).*z(1,1)]; % system of ode's

for i=1:3         % RK Startup
    k1 = f(t(i),y(:,i));
    k2 = f(t(i)+0.5*h,y(:,i)+0.5*h*k1);
    k3 = f((t(i)+0.5*h),(y(:,i)+0.5*h*k2));
    k4 = f((t(i)+h),(y(:,i)+k3*h));
    phi = (1/6)*(k1+2*k2+2*k3+k4);
     
    y(:,i+1) = y(:,i) + phi*h;
end
% vector to store function evaluations
fval = zeros(2,length(t));

fval(:,1) = f(t(1),y(:,1));
fval(:,2) = f(t(2),y(:,2));
fval(:,3) = f(t(3),y(:,3));
fval(:,4) = f(t(4),y(:,4));
% ABM method
for i = 4:length(t)-1;
    %Predictor
    y(:,i+1)= y(:,i) + ...
        (55 .* fval(:,i) - 59 .* fval(:,i-1) + ...
        37 .* fval(:,i-2) - 9 .* fval(:,i-3))...
        .* h24;
    %Evaluate
    fval(:,i+1) = f(t(i+1), y(:,i+1));
    %Corrector
    y(:,i+1) = y(:,i) + ...
        (9 .* fval(:,i+1) + 19 .* fval(:,i) - ...
        5 .* fval(:,i-1) + fval(:,i-2))...
        .* h24;
     %Evaluate
    fval(:,i+1) = f(t(i+1), y(:,i+1));
end

% now calculate the exact values
% Now to plot the true solution
% the true solution is a Bessel function:
%     y(t) = c1J1(t) + c2Y1(t)
%     y'(t) = c1(J0(t) - J2(t)) + c2(Y0(t) - Y2(t)),
%     with
%     d  = J1(1) [Y0(1) ? Y2(1)] ? Y1(1) [J0(1) ? J2(1)]
%     c1 = [Y0(1) ? Y2(1)] /d ? 1.36575994534763
%     c2 = ? [J0(1) ? J2(1)] /d ? ?0.51073987162512
t = 1:.1:2;
d = besselj(1,1).*(bessely(0,1)-bessely(2,1)) - ...
    bessely(1,1).*(besselj(0,1)-besselj(2,1));
c1 = (bessely(0,1) - bessely(2,1))./d;
c2 = -(besselj(0,1) - besselj(2,1))./d;
yb = c1*besselj(1,t) + c2*bessely(1,t);
pb = c1*(besselj(0,t) - besselj(2,t)) + c2*(bessely(0,t) - bessely(2,t));
figure
plot(t , yb , 'r*-' , t, pb , 'g*-' , t , y(1,:), 'mo-' , t , y(2,:) , 'co-')

% print out the results
disp(' '); disp('Solution Values'); disp(' '); disp(' t     y                y_exact          dydt        dydt_exact');...
    disp(' ================================================================');...
    for i = 1:length(t) disp(sprintf(' %4.2f%15.10f%15.10f%15.10f%15.10f',... 
        t(i),y(1,i),yb(i), y(2,i) , pb(i))); end
disp(' '); disp('Errors'); disp(' '); disp(' t       ABM           dydt_ABM');...
    disp(' =================================================');...
    for i = 1:length(t) disp(sprintf(' %4.2f%15.10f%15.10f%15.10f',t(i),y(1,i)-yb(i), y(2,i) - pb(i)... 
       )); end
