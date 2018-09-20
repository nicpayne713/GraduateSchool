x=[0:0.02:4];
y1=3*sin(pi.*x);
y2=exp(-.2.*x);
plot(x,y1,'r--')
hold on
plot(x,y2,'g--')
xlabel('x=0:02:4')
ylabel('3sin(pi*x),exp(-.2*x)')
title('Math 561 HW 1 Prob 4')
