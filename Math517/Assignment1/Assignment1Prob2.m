% Math 517
% Homework 1
% Question 1 AND 2
clear all
clc
% I have used a weighted Taylor expansion to sovle for u''(x2)
% Below is the system that has arisen
syms h1 h2 h3 
% O is the coefficient matrix for the 4x4 system of w0...w3
O = [1 1 1 1; 0 -h1 h2 h3+h2; 0 .5*h1^2 .5*h2^2 .5*(h2+h3)^2;
0 -(1/6)*h1^3 (1/6)*h2^3 (1/6)*(h2+h3)^3];
% Since we want to kill off everything except the u''(x2) term we use this
% v and solve O^-1 * v to sovle for the weights
v = [0 0 1 0]';

omega = inv(O)*v;

syms ux1 ux2 ux3 ux4 du4c1x du4c2x du4c3x du4um

error = omega(3)*(1/24)*h2^4*du4c2x + omega(4)*(1/24)*...
    (h2+h3)^4*du4c3x + omega(2)*(1/24)*h1^4*du4c1x;

error_collapsed = (1/24)*(omega(3)*h2^4+omega(4)*...
    (h2+h3)^4+omega(2)*h1^4)*du4um;

er = simplify(error_collapsed);

% Begin QUESTION 2
% function to compare error to
u = @sin;
x2 = pi/2;
% the true solution is u''(x2) = -sin(pi) = -1
true = -ones(length(H)); % vector for the loglog plot later on
H = linspace(1,1000,500);
du2dt_approx = [1, length(H)];

for i = 1:500;
        h = H(i)*[rand(1);rand(1);rand(1)]; %h(1) = h_1 etc.
        omega0 = -(2*(2*h(2) - h(1) + h(3)))/(h(1)*h(2)^2 + h(1)*h(3)*h(2));
        omega1 = (2*(2*h(2) + h(3)))/(h(1)^3 + ...
             2*h(1)^2*h(2) + h(3)*h(1)^2 + h(1)*h(2)^2 + h(3)*h(1)*h(2));
        omega2 =  (2*(h(2) - h(1) + h(3)))/(h(3)*h(2)^2 + h(1)*h(3)*h(2));
        omega3 = (2*(h(1) - h(2)))/(h(2)^2*h(3) + 2*h(2)*h(3)^2 +...
                h(1)*h(2)*h(3) + h(3)^3 + h(1)*h(3)^2);
            
        du2dt_approx(i) = omega0*u(x2) + omega1*u(x2-h(1)) + ...
                    omega2*u(x2+h(2)) + omega3*u(x2+h(2)+h(3));
end

% Plot the error against H on a log-log plot
figure(1)
loglog(H, du2dt_approx-true, 'r-o')
xlabel('H values for which 0<h1,h2,h3<H');
ylabel('Error of approximation vs. true solution')