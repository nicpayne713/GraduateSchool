%561 Final number 3
dydt = @(t,y) -10.*y+9.9.*exp(-1.*t./10);
%initial
y0 = 2;
h = .1; t = 0:h:2;
%Gear's third order method in the multivalue formulation
B = [1,1,1,1;0,1,2,3;0,0,1,3;0,0,0,1];
c = [6/11;1;6/11;1/11];
true = exp(-1.*t.*10) + exp(-1.*t./10);
dtrue = (-1/10).*exp(-1.*t.*10) + (-1/10).*exp(-1.*t./10);
d2true = (1/100).*exp(-1.*t.*10) + (1/100).* exp(-1.*t./10);
d3true = (-1/1000).*exp(-1.*t.*10) + (-1/1000).* exp(-1.*t./10);
y = zeros(4,length(t)); % initialize matrix to hold solutions
PEC  = [];
PEC2 = [];
PEC3 = [];
PEC4 = [];

% initial values
y(:,1) = [true(1);h.*dtrue(1);h.*d2true(1);h.*d3true(1)];

% PART a

%PEC
for i = 1:length(t)-1;
    % P
    y(:,i+1) = B*y(:,i);
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
end
PEC=y;
figure
plot(t,y(1,:),'r-',t,true,'o-')
title('Numerical and true solutions to ODE using PEC')
xlabel('t = 0:.1:2')
ylabel('outputs of numerical solution and true solution')
legend('numerical solution','true solution')

%P(EC)^2
for i = 1:length(t)-1;
    % P
    y(:,i+1) = B*y(:,i);
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;     
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
end
PEC2=y;
figure
plot(t,y(1,:),'g-',t,true,'o-')
title('Numerical and true solutions to ODE using P(EC)^2')
xlabel('t = 0:.1:2')
ylabel('outputs of numerical solution and true solution')
legend('numerical solution','true solution')

%P(EC)^3
for i = 1:length(t)-1;
    % P
    y(:,i+1) = B*y(:,i);
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
end

PEC3=y;
figure
plot(t,y(1,:),'b-',t,true,'o-')
title('Numerical and true solutions to ODE using P(EC)^3')
xlabel('t = 0:.1:2')
ylabel('outputs of numerical solution and true solution')
legend('numerical solution','true solution')

%P(EC)^4
for i = 1:length(t)-1;
    % P
    y(:,i+1) = B*y(:,i);
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
end
PEC4=y;
figure
plot(t,y(1,:),'c-',t,true,'o-')
title('Numerical and true solutions to ODE using P(EC)^4')
xlabel('t = 0:.1:2')
ylabel('outputs of numerical solution and true solution')
legend('numerical solution','true solution')

% all together
figure
plot(t, PEC(1,:),'r-')
hold on
plot(t,PEC2(1,:),'g-')
plot(t,PEC3(1,:),'b-')
plot(t,PEC4(1,:),'c-')
plot(t,true,'o')
hold off
title('solutions for P(EC)^m, m=1,2,3,4')
legend('PEC','P(EC)^2','P(EC)^3','P(EC)^4','true')
axis([0 2 .6 2])
% PART b

% If DG = (h*lambda,-1,0,0), and h*lambda = -1 then we have (from the
% stability notes online) that: 
W = 1./(1-(6/11).*(-1));

% PEC with W matrix
for i = 1:length(t)-1;
    % P
    y(:,i+1) = B*y(:,i);
    % EC
    y(:,i+1) = y(:,i+1) + W.*(h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
end
figure
plot(t,y(1,:),'r-',t,true,'o-',t,PEC3(1,:),'g',t,PEC4(1,:),'c')
title('Numerical and true solutions to ODE with W matrix in PEC')
xlabel('t = 0:.1:2')
ylabel('outputs of numerical solution and true solution')
legend('stabilized numerical solution','true solution', 'P(EC)^3','P(EC)^4')
