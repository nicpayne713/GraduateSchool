% b. classical 4th order Runge-Kutta
h=0.05;                      % step size
t = 1:h:1.7;                
y = zeros(1,length(t)); 
y(1) = 1;                     % initial condition
f_ty = @(t,y) y/t - t^2/y^2;  

for i=1:(length(t)-1)         % calculation loop
    k1 = f_ty(t(i),y(i));
    k2 = f_ty(t(i)+0.5*h,y(i)+0.5*h*k1);
    k3 = f_ty((t(i)+0.5*h),(y(i)+0.5*h*k2));
    k4 = f_ty((t(i)+h),(y(i)+k3*h));
    phi = (1/6)*(k1+2*k2+2*k3+k4);
    y(i+1) = y(i) + phi*h;  
end
y(length(t))
figure 
plot(t , y , '*', t, y, '-r')
title('4b Solution y(t) using RK4 method')
xlabel('Points where the slope of y is reevaluated')
ylabel('Outputs using RK4 method')
