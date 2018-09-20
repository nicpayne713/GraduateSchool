x = linspace(0,1,1000);
y = linspace(0,1,1000);
u = @(x,y) exp(x+0.5*y);
for p = 1:1000
    for q = 1:1000
        U(p,q) = u(x(p),y(q));
    end
end

figure;
mesh(x,y,U)

%%
clear all
close all
M517Assignment2Part1

figure;
mesh(x,y,usol);
title('Numerical solution');

figure;
mesh(x,y,Uhat);
title('True solution');

