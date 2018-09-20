% Math 517
% Homework 1
% Question 1
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