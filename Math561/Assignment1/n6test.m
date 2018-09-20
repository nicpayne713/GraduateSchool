upper =1 ;
lower = 0;
n=1; 
psum = 0; 
approxsum=0; 
while upper- lower > .0005 
    x =[1:n]; 
    psum = sum(((x.^(3/2) + x.^(1/2))).^(-1)); 
    upper = psum+ pi -2*atan(sqrt(n)) -1/(2*((n+1)^(3/2) + (n+1)^(1/2))); 
    lower = psum + pi -2*atan(sqrt(n+1)) +1/(2*((n+1)^(3/2) + (n+1)^(1/2)));
    n = n+1; 
end
upper
lower
n
