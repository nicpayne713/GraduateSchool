% Math 562 Assignment 1
% Problem 7 - FUNCTION

function [U ,S ,V] = my2dSVD(A)
% function for finding singular value decomposition of a 2x2 matrix, A
% Input: 2x2 matrix A
% Output: U - unitary matrix consisting of left singular vectors
%         V - unitary matrix consisting of right singular vectors
%         S - diagonal matrix with sigma_i along the diagonal
%
% function returns the 3 matrices and plots v_i, and u_i in the unit circle
% and ellipse as in Figure 4.1 on page 26 of Trefethen

global N Y c id s sigma1 sigma2 t u1 u2 v v1 v2

% sigma_1 = norm(A)_2 = max (over unit vectors) norm(Av)
N = exp(10); % number of partitions of unit circle
t = linspace(0, 2*pi, N);
% declare functions for finding singular vectors as plotting the unit
% circle
c = cos(t); s = sin(t);

% declare the vector to find sigma_1
v = [c;s]; % vector in polar coordinates
Y = A*v;
a = zeros(1,length(t));
for i = 1:N
    a(i) = norm(Y(:,i));
end
[sigma1, id] = max(a);
v1 = v(:,id);
u1 = (A*v1)/sigma1;
% now we rotate v1 and u1 by 90 degrees in order to extend each vector to
% an orthogonal basis for R^2
v2 = [0 -1;1 0]*v1;
u2 = [0 -1;1 0]*u1;
sigma2 = u2'*A*v2;
U = [u1 u2]
V = [v1 v2]
S = [sigma1 0; 0 sigma2]

% now for the plots

figure(1)

subplot(2,1,1)
plot(c,s)
hold on
arrow([0,0],v1,'b')
arrow([0,0],v2,'r')
title('Unit Circle with right singular vectors')
xlabel('x axis')
ylabel('y axis')
hold off

axis equal

subplot(2,1,2)
plot(Y(1,:),Y(2,:),'-')
hold on
arrow([0,0],u1,'b')
arrow([0,0],u2,'r')
title('Image of Unit Circle under A with left singular vectors')
xlabel('x axis')
ylabel('y axis')
hold off
end 
