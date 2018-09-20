function y = g(s)
%final homework number 4.
%y''(t)=3/2y^2(t), y(0)=4, y(1)=1
%PART a
% later s=linspace(-40,-5,100); 
%global y
y=zeros(1,length(s)); % initialize matrix to hold all values in
tspan=[0,1];
if nargin==0 %checks for argument
    warning('"g(s)" must contain an argument. Else, answer given on linspace.')
    return
end
format long
% ODE function
dy2dt = @(t,z) [z(2) ; (3/2).* z(1).^2];
for i = 1:length(s)
    [~ , Y] = ode45(dy2dt,tspan,[4 , s(i)]);

    y(i) = Y(end,1)-1; % find the error in the solution to y(1) 
    %based on initial guess y'(0) = s and the given initial value for y in the BVP
end
y' %returns the error for solution based on initial guess y'(0) = s, as a column vector
%varargout{y} = y;
end
