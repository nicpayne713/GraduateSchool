clear
clf
close all
q = 1;
w = @(x,y,t) (x-4.*t).^2 + (y-2.*t).^2;
x = linspace(-2,4);
y = linspace(-2,4);
[X,Y] = meshgrid(x,y);
R = zeros(length(x),length(y));
for t = 0:.005:.5
    L = w(X,Y,t);
    for i = 1:length(x)
        for j = 1:length(y)
            if L(i,j) <= 1
                R(i,j) = 1;
            end
        end
    end
    
	if t == 0
        figure(1)
        contour(X,Y,R)
        title('initial time')
    elseif t == .1
        figure(2)
        contour(X,Y,R)
        title('t1')
    elseif t == .25
        figure(3)
        contour(X,Y,R)
        title('t3')
    elseif t == .3
        figure(4)
        contour(X,Y,R)
        title('t4')
    elseif t == .4
        figure(5)
        contour(X,Y,R)
        title('t5')
        
    elseif t == .5
        figure(6)
        contour(X,Y,R)
        title('final time')
    end
end