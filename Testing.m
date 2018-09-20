% % Testing Document
% %%
% % AB 2-step
% f = @(r,z) r.^2 - (1+(3/2).*z).*r + (1/2).*z;
% theta = 3*pi/4;
% Rnew = 0;
% ds = .1;
% s =0;
% while abs(Rnew)<1
%     s = s+ds;
%     Rold = Rnew;
%     z = s.*exp(1i.*theta);
%     Rnew = max
%     
% end
% 
%     %% 
%     % Region of stability for 3rd order gear's method using PEC, P(EC)^2,
%     % etc.
%     %y' = -10y+9.9e^(-t/10)
%     % y(0) = 2
%     B = [1,1,1,1
%         0 1 2 3
%         0 0 1 3
%         0 0 0 1];
%     c = [6/11
%         1
%         6/11
%         1/11];
%     dydt = @(y,t) [y;
%         (.1)*(-10.*y+9.9.*exp(-t./10));
%         (.1^2)/2.*der((-10.*y+9.9.*exp(-t./10)));
%         (.1^3)/6.*der(der((-10.*y+9.9.*exp(-t./10))));];
%     t = 0:.1:2;
%     h=.1;
%     y = zeros(4,length(t));
%     y(:,1) = [2;0;0;0]; %initialize the y vector 
%     
%     %PEC
%     for i = 1:length(t);
%         % Predictor
%         y(:,i+1) = B .* y(:,i);
%         % Evaluate/Corrector
%         y(:,i+1) = y(:,i+1) + h.*(dydt(t(i+1),y(:,i+1)) - y(2,i+1)).*c;
%     end
%     
%     true = ode23s(dydt, t , 2);
%     plot(t , y(:,1) , 'ro-' , r , true , '*-')
%     legend('PEC' , 'True')
%     plot()
%     
%     
%%
% 561 Final homework #4
% BVP y''(t) = 3/2 y^2 (t) , y(0) = 4 , y(1) = 1
% part a

%%%%%%%%%%%%%%%%%%%%%% FUNCTION G IN MATH 561 FINAL FOLDER IS THE FUNCTION
%%%%%%%%%%%%%%%%%%%%%% THAT DOES WHAT THIS SCRIPT DOES!
dy2dt = @(t,z) [z(2) ; (3/2).* z(1).^2];
e_y = zeros(1,100);
tspan = [0,1];
g = @(s) ode45(dy2dt, tspan, [4 , s]);

% part b
s = linspace(-40,-5);
for i = 1:length(s)-1
    [T , Y] = g(s(i));
    e_y(2,i+1) = Y(length(Y'),1)-1; % find the error in the solution to y(1) based on initial guess y'(0) = s' and the given initial value for y in the BVP
end

plot( s , e_y , 'b', s, zeros(length(s)), 'c-') %y(1,:) is all zeros since nothing was caluclated here for y'
title('Error in solution to IVP for variable guess of solution to yprime(s) = 1')

indices = find([0 diff(sign(y(2,:)))]~=0)
% 2 , 14 , 93 are the indices where e_y changes sign - this is where the
% error is closest to 0.
% index 2 can be ignored since the graph starts at 0, so that counts as a
% sign change. But between s(13) and s(14) there is a sign change, and same
% for between s(92) and s(93). 
% So we'll use fzero to find the exact values of s to default accuracy

% can't use fzero on g since g doesn't actually find the difference. I
% either need to write g differently or use a for loop on values between
% s(13) and s(14) and another on on values between s(92) and s(93) to
% minimise the error then pick the best one possibly?

%%
%561 Final #5 GUIDE (This is only for part of the problem)
load sound.dat
z = sound(:,1);
Cz = sound(:,2);
cSpline = spline(z,Cz) % Spline the data to get a curve describing the
% speed of a sound wave based on the depth

coefs = cSpline.coefs;
cPcoefs = []; % initialze a matrix to hold coefficients of c prime

for p = 1:size(coefs,1) % this creats the coefficient matrix for the c prime spline
    X = coefs(p,:);
    cPcoefs = [ cPcoefs; polyder(X) ];
end
% first we mkpp (make a piecewise polynomial) out of the cPcoefs matrix
% using the mkpp. the breaks are the breaks of the cspline (data) and we
% use the matrix we just made.
cPrimeSpline = mkpp(cSpline.breaks, cPcoefs);
%now we make a function out of the spline to get information about any
%given depth z

c = @(z)ppval(cSpline,z);
cPrime = @(z)ppval(cPrimeSpline,z);

c(2000)
cPrime(2000)
%%
% 561 Final number 2
%b
syms h lambda
I = [1,0;0,1];
c = [1/2;2];
c2 = [0;1];
DG = [h*lambda,-1];
B = [1,1;0,1];

%PEC
Spec = (I + c*DG)*B

%PECE
Spece = (I + c2*DG)*(I + c*DG)*B

% check for zero stability by using h*lambda = 0
% We know from the notes that the region of stability extends for the
% entire negative real axis

subs(Spec,h*lambda,0)
subs(Spece,h*lambda,0)

% Spec has eigen values of 1 and -1, whose magnitidues are less or equal 1
% therefore it is zero stable, and Spece has eigen values 1 and 0, whose
% magnitudes are less than 1 therefore it's also zero stable.

%%
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
% matrix to store true values
TRUE = [true;dtrue;d2true;d3true];
y = zeros(4,length(t)); % initialize matrix to hold solutions
% initial values
y(:,1) = [true(1);h.*dtrue(1);h.*d2true(1);h.*d3true(1)];

%PEC
for i = 1:length(t)-1;
    % P
    y(:,i+1) = B*y(:,i);
    % EC
    y(:,i+1) = y(:,i+1) + (h.*dydt(t(i+1),y(1,i+1))-y(2,i+1)).*c;
end
figure
plot(t,y(1,:),'b-',t,true,'r-')
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
figure
plot(t,y(1,:),'b-',t,true,'r-')
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
figure
plot(t,y(1,:),'b-',t,true,'r-')
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
figure
plot(t,y(1,:),'b-',t,true,'r-')
title('Numerical and true solutions to ODE using P(EC)^4')
xlabel('t = 0:.1:2')
ylabel('outputs of numerical solution and true solution')
legend('numerical solution','true solution')
%Stability
%PEC
% Spec = (I + c*DG)*B
% 
% %PECEC
% Spec2 = (I + c*DG).^2*B
% 
% %P(EC)^3
% Spec3 = (I + c*DG).^3*B
% 
% %P(EC)^4
% Spec4 = (I + c*DG).^4*B

%%
% Possibly a code for number 1 on 561 FINAL HW
% Regions of stability

C = 'color'; c = {'b','r','g','m','y','c'};
x = [0 0]; y = [-8 8]; K = 'k'; LW = 'linewidth'; FS = 'fontsize';
plot(y,x,K,LW,1), hold on, plot(x,y,K)
t = chebfun(t,3,[0 2*pi]);
z = exp(1i*t); r = z-1;
s = 1; plot(r./s,C,c{1},LW,2)                      % order 1
s = (3-1./z)/2; plot(r./s,C,c{2},LW,2)             % order 2
s = (23-16./z+5./z.^2)/12; plot(r./s,C,c{3},LW,2)  % order 3
axis([-2.5 .5 -1.5 1.5]), axis square, grid on
title('Adams-Bashforth orders 1,2,3',FS,16)


%%

% Problem 1 on FINAL 561 hw
%global charpoly z
% the boundary for theta = pi/2 for the AB3 method is a little less than
% -.5 so I'll use that as my initial guess
% characteristic polynomial for AB4 is
charpoly = @(z) r.^4 - r.^3 - (55/14).*(z).*r.^2 - (37/24).*z.*r + 9.* z;
h = -pi/12;
for theta = pi/2:h:pi
end

    
%%
theta =linspace(0,2*pi,200);
z = exp(1i.*theta);
nu = (z.^4-z.^3)./((55/24)*z.^3 - (59/24)*z.^2 + (37/24)*z - 9/24);

plot(real(nu) , imag(nu)), axis equal, grid on

%%
theta = linspace(0,2*pi,200);
z=exp(1i.*theta);
hl=(24*(z.^4-z.^3))./(55*z.^3-59*z.^2+37*z-9);
plot(real(hl),imag(hl)),axis equal, grid on
%%
syms x
c = [1/2; 1];
c2 = [0 ;1];
dg = [x ,-1];
B = [1 1 ; 0 1];
I2 = eye(2);
S1 = (I2 + c*dg)*B;
S2 = (I2 + c2*dg)*(I2 + c*dg)*B;
lam1 = 0;
lam2 = 0;
z=[];
while isempty(find(abs(double(z))>1))
z = eig(subs(S1,lam1));
lam1 = lam1 - .005;
end
z=[];
while isempty(find(abs(double(z))>1))
z = eig(subs(S2,lam2));
lam2 = lam2 - .005;
end
lam1,lam2

%%
% 561 1c final homework
%1 C Brute Force
k = 4;
n=100;
AMcoef = [3/8,19/24,-5/24,1/24];
ABcoef = [ -3/8 , 37/24, -59/24 , 55/24];
theta = pi:-pi/(2*n):pi/2; 
tol = .5*exp(-12);
rr = zeros(n+1,1);
figure
hold on
for i = 1:n+1
    r=0;
    rho=0;
    while rho<1
        r = r+.01;
        z= r*exp(1i*theta(i));
        p = [ 1, -1-z*(AMcoef(1) + AMcoef(2)+z*AMcoef(1)*ABcoef(4)),...
    -z*(AMcoef(3:4)+z*AMcoef(1)*ABcoef(3:-1:2)), -z^2*AMcoef(1)*ABcoef(1)];
        rho = max(abs(roots(double(p))));
          
    end
    plot(r*cos(theta(i)),r*sin(theta(i)),'ro-',...
     r*cos(theta(i)),-r*sin(theta(i)),'ro-',...
     'LineWidth',2)
    grid on
    axis equal
end
%%
% Guide for Michael M for 561 #5 a
%561 Final #5 GUIDE (This is only for part of the problem)

z = [1,2,3,4,5,6,7,8];
Cz = [10,14,12,-10,2,8,20,-32];
cSpline = spline(z,Cz) % Spline the data to get a curve describing the
% relationship between z and Cz

coefs = cSpline.coefs;
cPcoefs = []; % initialze a matrix to hold coefficients of c prime

for p = 1:size(coefs,1) % this creats the coefficient matrix for the c prime spline
    X = coefs(p,:);
    cPcoefs = [ cPcoefs; polyder(X) ];
end
% first we mkpp (make a piecewise polynomial) out of the cPcoefs matrix
% using the mkpp. the breaks are the breaks of the cspline (data) and we
% use the matrix we just made.
cPrimeSpline = mkpp(cSpline.breaks, cPcoefs);
%now we make a function out of the spline to get information about any
%given depth z

c = @(z)ppval(cSpline,z);
cPrime = @(z)ppval(cPrimeSpline,z);

c(1) % should be 10 by construction of the problem
%%
% trying part 1 c of final 561 hw again

%global cpoly z
%z = exp(1i.*theta);
syms a b t
char = 2*a^3+(1/12)*(9*a*b*a^4 + 74*a*b*a^3-64*a*a^2*b + 38*a*b*a - 9*a*b);

cpoly = (19.*t.^2.*z)./6 - (3.*t.*z)./4 - (16.*t.^3.*z)./3 + (37.*t.^4.*z)./6 + (3.*t.^5.*z)./4 + 2.*t.^3;
n = 5;
for theta = pi:-pi/n:pi/2
    rr = rho(theta,-.5);
end
plot(rr.*cos(degtorad(theta)),rr.*sin(degtorad(theta)),'o')


%%
% Find the region of stability of a method numerically
% Follow AB-1 PEC through 4 steps
% region of stability changes every time
% it is not clear where this will settle down in the end

%global cpoly z
syms yn yn1 yn2 yn3 yn4 fn fn1 fn2 fn3 fn4 h lambda z r
n = 50;  % number of points
theta = pi:-pi/(2*n):pi/2;
theta(end+1) = 0;

% AB-1
fn = lambda * yn;
yn1 = yn + h*fn;         % P
fn1 = lambda * yn1;      % E
yn1 = yn + h/2*(fn1+fn); % C
%fn1 = lambda * yn1;      % E

% PEC Step 1
disp(' '); disp('PEC Step 1'); disp('===================');
E1 = subs(expand(yn1),h*lambda,z);
E1 = subs(E1,{yn},{1});
[cpoly, powers] = coeffs(r-E1,r)

rho_pec = zeros(size(theta));
rho0 = 3;
for j = 1:n+1
    rho_pec(j) = rho(theta(j),rho0);
    rho0 = rho_pec(j);
end
figure(1);
plot(rho_pec.*cos(theta),rho_pec.*sin(theta),'r',...
     rho_pec.*cos(theta),-rho_pec.*sin(theta),'r',...
     'LineWidth',2)
title('PEC Step 1','FontSize',20)
axis equal

% PEC Step 2
yn2 = yn1 + h*fn1;         % P
fn2 = lambda * yn2;        % E
yn2 = yn1 + h/2*(fn2+fn1); % C
%fn2 = lambda * yn2;        % E

disp(' '); disp('PEC Step 2'); disp('===================');
E2 = subs(expand(yn2),h*lambda,z);
E2 = subs(E2,{yn},{1});
[cpoly, powers] = coeffs(r-E2,r)

rho_pec2 = zeros(size(theta));
rho0 = 3;
for j = 1:n+1
    rho_pec2(j) = rho(theta(j),rho0);
    rho0 = rho_pec2(j);
end
figure(2);
plot(rho_pec2.*cos(theta),rho_pec2.*sin(theta),'g',...
     rho_pec2.*cos(theta),-rho_pec2.*sin(theta),'g',...
     'LineWidth',2)
title('PEC Step 2','FontSize',20)
axis equal

% PEC Step 3
yn3 = yn2 + h*fn2;         % P
fn3 = lambda * yn3;        % E
yn3 = yn2 + h/2*(fn3+fn2); % C
%fn3 = lambda * yn3;        % E

disp(' '); disp('PEC Step 3'); disp('===================');
E3 = subs(expand(yn3),h*lambda,z);
E3 = subs(E3,{yn},{1});
[cpoly, powers] = coeffs(r-E3,r)

rho_pec3 = zeros(size(theta));
rho0 = 3;
for j = 1:n+1
    rho_pec3(j) = rho(theta(j),rho0);
    rho0 = rho_pec3(j);
end
figure(3);
plot(rho_pec3.*cos(theta),rho_pec3.*sin(theta),'b',...
     rho_pec3.*cos(theta),-rho_pec3.*sin(theta),'b',...
     'LineWidth',2)
title('PEC Step 3','FontSize',20)
axis equal

% PEC Step 4
yn4 = yn3 + h*fn3;         % P
fn4 = lambda * yn4;        % E
yn4 = yn3 + h/2*(fn4+fn3); % C
%fn4 = lambda * yn4;        % E

disp(' '); disp('PEC Step 4'); disp('===================');
E4 = subs(expand(yn4),h*lambda,z);
E4 = subs(E4,{yn},{1});
[cpoly, powers] = coeffs(r-E4,r)

rho_pec4 = zeros(size(theta));
rho0 = 3;
for j = 1:n+1
    rho_pec4(j) = rho(theta(j),rho0);
    rho0 = rho_pec4(j);
end
figure(4);
plot(rho_pec4.*cos(theta),rho_pec4.*sin(theta),'k',...
     rho_pec4.*cos(theta),-rho_pec4.*sin(theta),'k',...
     'LineWidth',2)
title('PEC Step 4','FontSize',20)
axis equal

% all together
figure(5); clf;
plot(rho_pec.*cos(theta),rho_pec.*sin(theta),'r',...
     rho_pec.*cos(theta),-rho_pec.*sin(theta),'r',...
     rho_pec2.*cos(theta),rho_pec2.*sin(theta),'g',...
     rho_pec2.*cos(theta),-rho_pec2.*sin(theta),'g',...
     rho_pec3.*cos(theta),rho_pec3.*sin(theta),'b',...
     rho_pec3.*cos(theta),-rho_pec3.*sin(theta),'b',...
     rho_pec4.*cos(theta),rho_pec4.*sin(theta),'k',...
     rho_pec4.*cos(theta),-rho_pec4.*sin(theta),'k')
%      rho_pec_mv.*cos(theta),rho_pec_mv.*sin(theta),'m',... % limiting region
%      rho_pec_mv.*cos(theta),-rho_pec_mv.*sin(theta),'m',...% from multivalue method
     %'LineWidth',2);
title('PEC, Steps 1 - 4','FontSize',20)
axis equal
%%
% Find the region of stability of a method numerically
% This is for ABM 4th order methodPECE

global cpoly z
syms yn yn1 yn2 yn3 yn4 fn fn1 fn2 fn3 fn4 h lambda z r
n = 50;  % number of points
theta = pi:-pi/(2*n):pi/2;
theta(end+1) = 0;

% % ABM 4-step method in multivalue form
% fn = lambda * yn;
% yn1 = yn + h*fn;         % P
% fn1 = lambda * yn1;      % E
% yn1 = yn + h/2*(fn1+fn); % C
% fn1 = lambda * yn1;      % E

B = [1,55/24,-59/24,37/24,-9/24;  % AB 4-step predictor
     0,4,-6,-4,1;
     0,1,0,0,0;
     0,0,1,0,0;
     0,0,0,1,0];  
c = [9/24,];  % AB 4th order corrector
c2 = [0; 1; ; ;0];  % final evaluation
DG = [z,-1,0,0];
I = eye(5);

% PEC
S = (I + c*DG)*B
cpoly = charpoly(S);

rho_pec = zeros(size(theta));
rho0 = 3;
for j = 1:n+1
    rho_pec(j) = rho(theta(j),rho0);
    rho0 = rho_pec(j);
end
figure(1);
plot(rho_pec.*cos(theta),rho_pec.*sin(theta),'r',...
     rho_pec.*cos(theta),-rho_pec.*sin(theta),'r',...
     'LineWidth',2)
title('PEC','FontSize',20)
axis equal

% PECE
S = (I + c2*DG)*(I + c*DG)*B
cpoly = charpoly(S);

rho_pece = zeros(size(theta));
rho0 = 3;
for j = 1:n+1
    rho_pece(j) = rho(theta(j),rho0);
    rho0 = rho_pece(j);
end
figure(2);
plot(rho_pece.*cos(theta),rho_pece.*sin(theta),'g',...
     rho_pece.*cos(theta),-rho_pece.*sin(theta),'g',...
     'LineWidth',2)
title('PECE','FontSize',20)
axis equal

% all together
figure(3); clf;
plot(rho_pec.*cos(theta),rho_pec.*sin(theta),'r',...
     rho_pec.*cos(theta),-rho_pec.*sin(theta),'r',...
     rho_pece.*cos(theta),rho_pece.*sin(theta),'g',...
     rho_pece.*cos(theta),-rho_pece.*sin(theta),'g',...
     'LineWidth',2);
title('PEC and PECE','FontSize',20)
axis equal
%%
g = 9.81; m = 68.1; t = 12; cd = 0.25;
v = sqrt(g * m / cd) * tanh(sqrt(g * cd / m) * t)

%%
% Math 517
% Homework 1
% Questions 1-3
clear all
clc
% I have used a weighted Taylor expansion to sovle for u''(x2)
% Below is the system that has arisen
syms h1 h2 h3 
% O is the coefficient matrix for the 4x4 system of w0...w3
O = [1 1 1 1; 0 -h1 h2 h3+h2; 0 .5*h1^2 .5*h2^2 .5*(h2+h3)^2;
0 -(1/6)*h1^3 (1/6)*h2^3 (1/6)*(h2+h3)^3];
% Since we want to kill off everything except the u''(x2) term we use this
% v and solve O^-1 * v to sovle for the weights
v = [0 0 1 0]';

omega = inv(O)*v;
omega = simplify(omega)

syms ux1 ux2 ux3 ux4 du4c1x du4c2x du4c3x du4um

error = omega(3)*(1/24)*h2^4*du4c2x + omega(4)*(1/24)*...
    (h2+h3)^4*du4c3x + omega(2)*(1/24)*h1^4*du4c1x;

error_collapsed = (1/24)*(omega(3)*h2^4+omega(4)*...
    (h2+h3)^4+omega(2)*h1^4)*du4um;

er = simplify(error_collapsed);

% Begin QUESTION 2
% function to compare error to
u = @sin;
x2 = pi/2;
H = linspace(exp(-4),exp(-2),500);
% the true solution is u''(x2) = -sin(pi) = -1
true = -ones(1,length(H)); % vector for the loglog plot later on
du2dt_approx = [1, length(H)];

for i = 1:500;
        h = H(i)*rand(1); %h(1) = h_1 etc.
        omega0 = -(2*(2*h - h + h))/(h*h^2 + h*h*h);
        omega1 = (2*(2*h + h))/(h^3 + ...
             2*h^2*h + h*h^2 + h*h^2 + h*h*h);
        omega2 =  (2*(h - h + h))/(h*h^2 + h*h*h);
        omega3 = (2*(h - h))/(h^2*h + 2*h*h^2 +...
                h*h*h + h^3 + h*h^2);
            
        du2dt_approx(i) = omega0*u(x2) + omega1*u(x2-h) + ...
                    omega2*u(x2+h) + omega3*u(x2+h+h);
end

% Plot the error against H on a log-log plot
errorHi = abs(du2dt_approx - true);
figure(1)
loglog(H, errorHi, 'bo')
%set ( gca, 'xdir', 'reverse' )

xlabel('H values for which 0<h1,h2,h3<H');
ylabel('Error of approximation vs. true solution')

% BEGIN QUESTION 3
% Least-squares fit: log(E(H)) = K + plog(H)

L = [(log(du2dt_approx))]';
A = [ ones(500,1) , L ];
b = [(log(errorHi))]';
% x = [K p]
x = (A'*A)\A'*b
abs(x)
% y = pinv(A)*b
% abs(y)

%%
% 533 HW 1 #2
%A
clear all
clc
txt = 'MSCUPTEHCIPHCJICPTEHCQBPUSPJKEOEHCKSCRNCQVOEHCQYXKDCKTEHCJLESXQQVLCEQLCPUGFEICSEICFEICIEQCQPONGDVKSPUMNDJEUSPQGNULFEEDNGLCSEEBCQQEEJCBCDKKQXIPUCNJCUSEMCQHCUUCGVTESUJEUUNWEUREQGFPQKCQQCPLPUYXJSFDNUKEUYXJSEYPIP';
ychar=zeros(26,1);
yct=zeros(26,1);
disp('Frequency of each character')
for k=0:25,
   ychar(k+1)=char('A'+k);
   v=(txt==ychar(k+1,1));   %perform check to see 
   yct(k+1)=sum(v);
   ydisp=[ychar(k+1), '  ',num2str((yct(k+1))/length(txt))];
   disp(ydisp)
end
%end;
IOC = sum((yct./length(txt)).^2)*26; % Calculcates the Index of Coincidence
disp('IOC = '), disp(IOC) % = 
Incidence = 0;
for i = 0:25
    p_i = yct(i+1)/length(txt);
    p_i = (p_i)^2;
    Incidence = p_i + Incidence;
end
disp(Incidence*26) % I keep getting 4.4--- WHY?

%B
% Look for common digrams in the ciphertext
% there are 26^2 possible digrams
digramtab = zeros(26,26);
% for k = 1:26
%     for l = 1:26
%         if 



%%
% Math 517 Homeowrk 5 Problem 5
% Applying methods to given PDE
% 
%
% First set of initial conditions: u(0,x) = 2exp(-200(x-.5)^2)
% Second set: u(0,x) = 2exp(-200(x-.5)^2)cos(40pi*x)

h = [1/10,1/20,1/40, 1/80]; % spacial step size
Tfinal = 1; % final time
a = 1;

% Initialize vector to hold error
ErrorLW = zeros(1,length(h));
ErrorRK = zeros(1,length(h));
AccuracyLW = zeros(1,length(h));
AccuracyRK = zeros(1,length(h));

for i = 1:length(h) % big loop for accuracy analysis
    M = 1/h(i) - 1; % number of intervals
    k = h(i); % time step
    x = linspace(0,1,M+2); % grid points
    T = 0:k:Tfinal; % time grid points
    
    % ut = -aux advection equation
    U = zeros(length(T),length(x)); % intialize U matrix
    Ur = zeros(length(T),length(x)); % initialize Ur matrix for RK3 method
    
    % initial conditions
    % 1)
%     U(1,:) = 2.*exp(-200.*(x-.5).^2);
%     Ur(1,:) = 2.*exp(-200.*(x-.5).^2);
%     Exact = 2*exp(-200*(x(end)-T(end)-.5)^2); % exact solution at final time
    % 2)
    U(1,:) = 2.*exp(-200.*(x-.5).^2).*cos(40*pi.*x);
    Ur(1,:) = 2.*exp(-200.*(x-.5).^2).*cos(40*pi.*x);
    Exact = 2.*exp(-200.*(x-Tfinal-.5).^2).*cos(40*pi.*(x-Tfinal));
    
    % Method from Question 3 (LxW type)
    
    % Big time loop
    for t = 2:length(T)     
        for j = 1:M+2 % spacial loop
            jm2 = j-2; jm1 = j-1; jp1 = j+1;
            if jm2 < 1
                jm2 = jm2+M+2;
            end
            if jm1 < 1
                jm1 = jm1+M+2;
            end
            if jp1 > M+2
                jp1 = jp1-M-2;
            end
            U(t,j) = U(t-1,j) - (k*a)/(6*h(i))*(U(t-1,jm2)-6*U(t-1,jm1)+3*U(t-1,j)+2*U(t-1,jp1)) + ...
                (k^2*a^2)/(2*h(i)^2)*(U(t-1,jm1)-2*U(t-1,j)+U(t-1,jp1)) - ...
                (k^3*a^3)/(6*h(i)^3)*(-U(t-1,jm2)+3*U(t-1,jm1)-3*U(t-1,j)+U(t-1,jp1));
        end 
    end
    
    % Method from Question 4 (RK3 type)
    
    % Big time loop
    for t = 2:length(T)
        % Loop for Y1
        Y1 = zeros(1,length(x));
        Us = zeros(1,length(x)); % Ustar vector
        for j = 1:M+2 % spacial loop
            jm2 = j-2; jm1 = j-1; jp1 = j+1;
            if jm2 < 1
                jm2 = jm2+M+2;
            end
            if jm1 < 1
                jm1 = jm1+M+2;
            end
            if jp1 > M+2
                jp1 = jp1-M-2;
            end
            Y1(j) = -a/(6*h(i))*(Ur(t-1,jm2) - 6*Ur(t-1,jm1)+3*Ur(t-1,j)+2*Ur(t-1,jp1));
            Us(j) = Ur(t-1,j)+.5*k*Y1(j);
        end
        
        % Loop for Y2
        Y2 = zeros(1,length(x));
        Uss = zeros(1,length(x)); % Ustarstar vector
        for j = 1:M+2 % spacial loop
            jm2 = j-2; jm1 = j-1; jp1 = j+1;
            if jm2 < 1
                jm2 = jm2+M+2;
            end
            if jm1 < 1
                jm1 = jm1+M+2;
            end
            if jp1 > M+2
                jp1 = jp1-M-2;
            end
            Y2(j) = -a/(6*h(i))*(Us(jm2) - 6*Us(jm1)+3*Us(j)+2*Us(jp1));
            Uss(j) = Ur(t-1,j)+.75*k*Y2(j);
        end
        
        % Loop for Y3
        Y3 = zeros(1,length(x));
        for j = 1:M+2 % spacial loop
            jm2 = j-2; jm1 = j-1; jp1 = j+1;
            if jm2 < 1
                jm2 = jm2+M+2;
            end
            if jm1 < 1
                jm1 = jm1+M+2;
            end
            if jp1 > M+2
                jp1 = jp1-M-2;
            end
            Y3(j) = -a/(6*h(i))*(Uss(jm2) - 6*Uss(jm1)+3*Uss(j)+2*Uss(jp1));
        end
        
        % Update Loop
        for j = 1:M+2
            Ur(t,j) = Ur(t-1,j)+(k/9)*(2*Y1(j) + 3*Y2(j) + 4*Y3(j));
        end
    end
    
    % Error analysis 
    ErrorLW(i) = norm(U(M+2,:) - Exact,2)/norm(Exact,2);
    ErrorRK(i) = norm(Ur(M+2,:) - Exact,2)/norm(Exact,2);
    
end

% Accuracy Analysis (use log 2 since h is halved in each successive step)
for i = 2:length(h)
    AccuracyLW(i) = log2(ErrorLW(i-1)/ErrorLW(i));
    AccuracyRK(i) = log2(ErrorRK(i-1)/ErrorRK(i));
end

Ta = table(h', AccuracyLW', AccuracyRK', 'VariableNames',{'stepsize' 'LWAccuracy' 'RKAccuracy'});
disp(Ta)

[XX,TT] = meshgrid(x,T);
Ex = 2.*exp(-200.*(XX-TT-.5).^2).*cos(40*pi.*(XX-TT));

figure
mesh(x,T,U)
title('Numerical LW')
xlabel('x')

figure
mesh(x,T,Ex)
title('Exact')
xlabel('x')

figure
mesh(x,T,Ur)
title('Numerical RK')
xlabel('x')



































