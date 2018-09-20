% Problem 4
% Part C
% Extra Credit
%
t = linspace(0,2*pi);
syms z
a = ones(100,1);
b = (2*z+6)/(z^2-4*z+6) * a;
c = b - exp(1i*t)';
for k = 1:100;
    plot(solve(c(k)),'o')
    hold on
end
title('Region of Stability for given implicit 2 stage-RK method')