
function F = mjr(z)
a = z(1);
b = z(2);
c = z(3);
d = z(4);
F(1) = (4.608*10^(6)*c)/b - 800 ;
F(2) = 12.5*((4.608*10^(6)*(d-c)/b))-16000;
F(3) = (10/3)*c+12.5*a*(d-c)^2-b;
F(4) = c*(10*c+12.5*a)-.5*c*(10-c)-(d-c)*12.5*a;
end







