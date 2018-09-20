%Nicholas Payne
%Math 561 Homework 1
%Due September 11, 2014

close all
clear all
clc
format long
%%
%---------------------------Problem 2-------------------------------
%The largest factorial that can be computed is the 1 which is less or equal
%to the solution of equation (1.6) on page 4 of Gautschi for t=23, s=7, which is 2^127. 
%It turns out that 34! is greater than this and 33! is less, so 33 is our
%upperbound for the loop below.
maxfa=1;
for i=1:33 %Upperbound
    if factorial(i) <= 2^(52) %Largest value able to be converted from base 
        %10 to base 2
        r = dec2base(factorial(i),2);
        
    else
        i= 33; %This will end the loop 
    end

    %Since matlab is double-precision I need to make sure that the number I
    %got is actually representable in single precision
    %Need to make sure that there are no 1's after the 23rd place in the
    %Mantissa

for k=17:-1:1 % 17 is the largest factorial less than 2^52, found in the loop above
    n = find(num2str(r)=='1'); %Finds the 1's in the binary rep of the factorial
    if ~isempty(n)
        if (n(end)-n(1)) <= 23 %Checks that there are no 1's beyond the 23rd bit
            
        else %If there are then the program incriments down and checks the next natural number
            r=dec2base(factorial(k-1),2);
            maxfa=k-1;
        end
    end
end
end
MaxFactorialInIEEESinglePrecision=maxfa;
MaxFactorialInIEEESinglePrecision

%%
%---------------------------Problem 3-------------------------------
%PART A

p=[5,3,-1,0]; q=[1,-1,2];

%multiplying p and q is given by
prod = conv(p,q);
product = poly2sym(prod, 'x');
product

%to divide the poly nomials use the deconv() command
[Q,R] = deconv(p,q);
quotient = poly2sym(Q,'x');
remainder = poly2sym(R,'x');
quotient
remainder


%to find the derivative of a polynomial, use polyder(p)
der = polyder(p);
derivative = poly2sym(der, 'x');
derivative

%PART B

%polyadd written as a function in separate script file
s = polyadd(p,q); %prints the answer as a polynomial 

strcat('Adding the polynomials yields: ',char(s)) 

%PART C

%Writing polyadd as an inline function
f = inline('[zeros(1,m-length(f)) f] + [zeros(1,m-length(g)) g]','f','g','m'); 
f(p,q,5) %adds the 2 polynomials using the inline function

%%
%-------------------------Problem 4----------------------------------------
x=[0:0.02:4];
y1=3*sin(pi.*x);
y2=exp(-.2.*x);
plot(x,y1,'r--')
hold on
plot(x,y2,'g--')
xlabel('x=0:02:4')
ylabel('3sin(pi*x),exp(-.2*x)')
title('Math 561 HW 1 Prob 4')
gtext('Intersection')
%Plot included after MATLAB Outputs section
%%
%-------------------------Problem 5----------------------------------------
%PARTS A and C are included after the MATLAB Outputs section

%PART B
%I defined a function in another script as follows
%function y = m561_1_5(a,b) 
%for i=1:99
%y = ((5/2)*b -a); 
%a=b;
%b=y;
%end
%end

%For our first set of initial values y_0 = 1, y_1=1/5
FirstInitialValuesy_100 = m561_1_5(1,.5);
FirstInitialValuesy_100

%For our second set of initial values y_0=1/3, y_1=1/6
SecondInitialValuesy_100 = m561_1_5(1/3,1/6);
SecondInitialValuesy_100


%%
%----------------------------Problem 6 (EC)--------------------------------
ups = 1; %k=1 uppersum from bound
lws = 0; %k=0 lowersum from bound
n=1;
pasum = 0; %partial sum
asum=0;
while ups- lws > .0001; %safe bound of 4 decimal places
    x=[1:n];
    pasum = sum(((x.^(3/2) + x.^(1/2))).^(-1)); 
    ups = pasum + pi - 2*atan(sqrt(n)) - 1/(2*((n+1)^(3/2)+(n+1)^(1/2)));
    lws = pasum + pi - 2*atan(sqrt(n+1)) + 1/(2*((n+1)^(3/2)+(n+1)^(1/2)));
    n=n+1;
end
ups
lws
n
