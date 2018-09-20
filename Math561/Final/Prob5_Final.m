% 561
% Final: Problem 5
clear all
close all
clc
format long
% PART a
% Need to be able to 'trace a ray' ie. 'solve the ODE and plot the
% resulting curve'. Start by tracing the ray beginning at z0 = 2000 and
% theta0 = 5.4 degrees from x = linspace[0,24] miles which is 
x =linspace(0,145824); % feet
% using the function raytrace (code included), we can solve g(x) where x is the
% linspace above
% spline the sound.dat
load sound.dat
global A coefs cPcoefs cSpline cPrimeSpline c cPrime dz2dx
z = sound(:,1);
Cz = sound(:,2);
cSpline = spline(z,Cz); % Spline the data to get a curve describing the
% speed of a sound wave based on the depth

coefs = cSpline.coefs;
cPcoefs = zeros(16,3); % initialze a matrix to hold coefficients of c prime

for p = 1:size(coefs,1) % this creats the coefficient matrix for the c prime spline
    X = coefs(p,:);
    cPcoefs(p,:) = polyder(X);
end
% first we mkpp (make a piecewise polynomial) out of the cPcoefs matrix
% using the mkpp. the breaks are the breaks of the cspline (data) and we
% use the matrix we just made.
cPrimeSpline = mkpp(cSpline.breaks, cPcoefs);
%now we make a function out of the spline to get information about any
%given depth z
c = @(z) ppval(cSpline,z);
cPrime = @(z) ppval(cPrimeSpline,z);

% Snell's Law
A = @(z,theta) cos(degtorad(theta))/c(z); 

%ODE
% y1 = z , y2 = z'
% y1' = y2 , y2' = -c'(y1)/(A^2*c(y1)^3)

dz2dx = @(t,y) [y(2);-cPrime(y(1))./((A(y(1),radtodeg(atan(y(2))))).^2.*c(y(1)).^3)];
% initial conditions
y = zeros(1,2);
y(1) = 2000;
y(2) = tan(degtorad(5.4));
[T,Y] = ode45(dz2dx,x, y);
figure
plot(T , Y(:,1) , 'b-', T , Y(:,2) , 'g-')
title('path or ray starting at 5.4 degrees')
xlabel('distance x[ft]')
ylabel('depth z[ft]')
axis ij



% PART b
% raytrace only takes 1 input, so this for loop calculates by how much the
% ray misses the receiver given each theta0
theta = -10:1:10;
error = zeros(1,21);
for i = 1:21
error(1,i) = raytrace(theta(i)); %raytrace is the function which measures
    % the vertical distance by which the ray starting with angle theta0
    % misses the receiver
end
Table = [theta',error'] %This prints out an array expressing the error for each theta0
%Table =
%    1.0e+03 *
% 
%   -0.010000000000000   1.864011296023142
%   -0.009000000000000   0.434600900416246
%   -0.008000000000000  -0.169628760986016
%   -0.007000000000000  -1.539259656434191
%   -0.006000000000000  -1.325382475472145
%   -0.005000000000000  -1.400530557160409
%   -0.004000000000000   0.448663998304038
%   -0.003000000000000   0.141263997722316
%   -0.002000000000000   0.172468769801948
%   -0.001000000000000   0.089980855906799
%                    0   0.258765951215780
%    0.001000000000000   0.410094789017413
%    0.002000000000000   0.662950242052646
%    0.003000000000000   0.872191065183585
%    0.004000000000000  -0.922391698865355
%    0.005000000000000  -0.466664006062790
%    0.006000000000000   0.513957346625004
%    0.007000000000000   0.387231061578991
%    0.008000000000000  -1.623178681265491
%    0.009000000000000  -1.847430848265012
%    0.010000000000000  -1.252962437485318

%PART c
theta1 = fzero(@raytrace,-8)
theta2 = fzero(@raytrace,7)
theta3 = fzero(@raytrace,4)

% theta1 =
% 
%   -8.371979844531827
% 
% 
% theta2 =
% 
%    7.240417364520359
% 
% 
% theta3 =
% 
%    3.760096932242514

w = zeros(1,2); w(1) = 2000 ; w(2);% = tan(degtorad(thetaN))
[T,Y1] = ode45(dz2dx,x,[w(1),tan(degtorad(theta1))]);
[~,Y2] = ode45(dz2dx,x,[w(1),tan(degtorad(theta2))]);
[~,Y3] = ode45(dz2dx,x,[w(1),tan(degtorad(theta3))]);
figure
plot(T , Y1(:,1) , 'b-', T , Y2(:,1) , 'g-' , T,Y3(:,1),'r-')
title('path of rays with optimal initial angles to hit the receiver')
xlabel('distance x[ft]')
ylabel('depth z[ft]')
legend('theta near -8 degrees','theta near 7 degrees','theta near 4 degrees')
axis ij

%PART d
% find how long it takes for the ray with initial angle theta3 to reach the
% receiver
% speed = distance/time so time = distance/speed so we integrate
% T = int[(sqrt(1+(z'(x))^2]/[c(z(x))] dx
% First we make functions z(x) and z'(x) in a similar way to how we made
% polynomials for c(z and c'(z) in part a.
% theta here is theta3 from part c
zSpline = spline (x,Y3(:,1));
zPrimeSpline = spline(x,Y3(:,2));
z = @(x) ppval(zSpline,x);
zPrime = @(x) ppval(zPrimeSpline,x);
h = @(x) (1+(zPrime(x)).^2).^(1/2)./(c(z(x)));
Time = quad(h,x(1),x(end))
% Time =
% 
%   29.951114880691527
