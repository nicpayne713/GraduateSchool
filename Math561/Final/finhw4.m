function g=finhw4(s);
%final homework number 4.
%y''(t)=3/2y^2(t), y(0)=4, y(1)=1
%This function can take an argument which is a 1xn or nx1 array.
%Part A
%s=linspace(-40,-5,100); 
g=[0];
tspan=[0,1];
if nargin==0 %checks for argument
    warning('Must have an argument in function "finhw4a()" Else, answer given on linspace.')
    return
end
format long
length(s);
    function dy=finprob4(t,y) %ODE function 
        dy=zeros(2,1);
        dy(1)=y(2);
        dy(2)= 3/2.*y(1).^2;
%     count_call = count_call + 1;
    end;
for n=1:length(s) %This function allows for input of a vector.
b=s(n);
x=zeros(2,1); x(1)=4; x(2)=b; %initial conditions
[t,y]=ode45(@finprob4,tspan,x); %function evaulation
%     hold on
end
%str=sprintf('y(t) given yprime(0)=%6.15f',s); %put full double pression variable in title
plot(t,y(:,1),'-r'); xlabel('t'); ylabel('y'); title(str); %graph y(t) given s
g(n)=y(end,1)-1; %returns g(s)
end