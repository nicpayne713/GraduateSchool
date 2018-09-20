% Assignment 5
% Problem 3
% PART b
%
clear all
close all
clc
h = 0.1; % stepsize
h24 = h/24;
t = 1:h:2; % interval
y = zeros(1, length(t)); % vector to store y values in
p = zeros(1, length(t)); % vector to store y' values in
f1 = @(t , y , p) p; % function for substitution
f2 = @(t , y , p) -p./t - (t.^2-1)./(t.^2).*y;
y(1) = 1; % initial conditions
p(1) = 0;
% begin with the RK4 method to start the ABM method
for i = 1:3;
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
% Now implement the ABM Method
for i = 4:length(t)-1;
    % Predictor
    y(i+1)= y(i) + ...
        (55 .* f1(t(i) , y(i) , p(i)) - 59 .* f1(t(i), y(i-1) , p(i-1)) + ...
        37 .* f1(t(i) , y(i-2) , p(i-2)) - 9 .* f1(t(i) , y(i-3) , p(i-3)))...
        .* h24;
    
    p(i+1)= p(i) + ...
        (55 .* f2(t(i) , y(i) , p(i)) - 59 .* f2(t(i), y(i-1) , p(i-1)) + ...
        37 .* f2(t(i) , y(i-2) , p(i-2)) - 9 .* f2(t(i) , y(i-3) , p(i-3)))...
        .* h24;
%     % Evaluate 
%     y(i+1) = f1(t(i+1), y(i+1), p(i+1));
%     p(i+1) = f2(t(i+1), y(i+1), p(i+1));
    % Corrector
    y(i+1) = y(i) + ...
        (9 .* f1(t(i+1), y(i+1) , p(i+1)) + 19 .* f1(t(i+1), y(i) , p(i)) - ...
        5 .* f1(t(i+1), y(i-1) , p(i-1)) + f1(t(i+1), y(i-2) , p(i-2)))...
        .* h24;

    p(i+1) = p(i) + ...
        (9 .* f2(t(i+1), y(i+1) , p(i+1)) + 19 .* f2(t(i+1), y(i) , p(i)) - ...
        5 .* f2(t(i+1), y(i-1) , p(i-1)) + f2(t(i+1), y(i-2) , p(i-2)))...
        .* h24;
%         % Evaluate 
%     y(i+1) = f1(t(i+1), y(i+1), p(i+1));
%     p(i+1) = f2(t(i+1), y(i+1), p(i+1));
end

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
yb = c1*besselj(1,t) + c2*bessely(1,t)
pb = c1*(besselj(0,t) - besselj(2,t)) + c2*(bessely(0,t) - bessely(2,t))
figure
plot(t , yb , 'r*-' , t, pb , 'g*-' , t , y, 'mo-' , t , p , 'co-')
% print out the results
y
p
exact_y = yb
exact_y_prime = pb
diff_y = y-yb
diff_y_prime = p-pb