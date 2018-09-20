function [D] = ode_euler(phi, h, y_0 , a , b)
y = y_0;
t = a+h;
while t <=b;
    y = y + h * phi;
    t = t+h;
    D = y;
end

    
