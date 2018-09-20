%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PART a)
%
clc
format long
f5 = @(x) exp(2.*x)./((1+x.^2).^2);
[q,n] = quad(f5, 0,3)
[ql,nl] = quadl(f5,0,3)
qgk = quadgk(f5,0,3, 'MaxIntervalCount', 670)
int = integral(f5, 0,3)
%
% PART b)
% In Tex File
% PART c)
%
fc = @(x,y) 2.*y.*sin(x)+(cos(x)).^2;
c = @(x) sin(x);
d = @(x) cos(x);
quad2d(fc,0,pi/4,c,d)
integral2(fc,0,pi/4,c,d)