% Final
% Problem 1
%%
% part a
clc
close all
clear all
% plot the region of stability for AB4
% find the characteristic equation for AB4, solve for hlambda, and let
% theta run from 0 to 2pi, plot the resulting points
theta =linspace(0,2*pi,500);
z = exp(1i.*theta); 
% nu is hlambda
nu = (z.^4-z.^3)./((55/24)*z.^3 - (59/24)*z.^2 + (37/24)*z - 9/24);
plot(real(nu) , imag(nu))
axis equal
hold on

% PART b

% ADAMS MOULTON 3
xi = (z.^3-z.^2).*24./(9*z.^3+19*z.^2-5.*z+1);
plot(real(xi),imag(xi),'r-')

%PART c
%
k = 4;
n=100;
AMcoef = [3/8,19/24,-5/24,1/24];
ABcoef = [ -3/8 , 37/24, -59/24 , 55/24];
theta = pi:-pi/(2*n):pi/2; 
tol = .5*exp(-14);
rr = zeros(n+1,1);

for i = 1:n+1
    r=0;
    rho=0;
    while rho<1
        r = r+.01;
        z= r*exp(1i*theta(i));
        p = [ 1, -1-z*(AMcoef(1) + AMcoef(2)+z*AMcoef(1)*ABcoef(4)),...
    -z*(AMcoef(3:4)+z*AMcoef(1)*ABcoef(3:-1:2)), -z^2*AMcoef(1)*ABcoef(1)];
        rho = max(abs(roots(double(p))));
          
    end
    plot(r*cos(theta(i)),r*sin(theta(i)),'rp',...
     r*cos(theta(i)),-r*sin(theta(i)),'rp',...
     'LineWidth',2)
    grid on
    axis equal
    axis([-3.5 2 -2 2])
    legend('AB4','AM3','4PECE')
    xlabel('Re')
    ylabel('Im')
end
%%
% part c a different way
% THIS DID NOT WORK - CHARACTERISTIC EQUATION IS WRONG
%alpha = -(24.*z.^3)./(38.*z - 9 - 64.*z.^2 + 74.*z.^3 + 9.*z.^4);
%alpha = (16.*z.^3)./(3 - 12.*z + 18.*z.^2 - 12.*z.^3 + 3.*z.^4);
%alpha = -(48.*z.^3)./(38.*z - 9 - 64.*z.^2 + 74.*z.^3 + 9.*z.^4);
%plot(real(alpha),imag(alpha),'g-')