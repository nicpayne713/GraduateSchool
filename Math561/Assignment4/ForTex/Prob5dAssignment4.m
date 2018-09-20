% Problem 5d
% Bessel's Equation of Order One
%
% declare an array for the MATLAB besselj function to evaluate 
% the ODE at such that it's evaluated at the same points as in Euler's
% Method and with the RK4 Method
z = [1:.1:2] ;
J = besselj(1, z)
figure
plot(z,J,'o', z, J, 'r')
title('Graph of Bessels Equation of the first kind of order one')
xlabel('Points used for solving the ODE with Euler and RK4 using h = .1')
ylabel('Outputs according the the MATLAB bessel equation solver')