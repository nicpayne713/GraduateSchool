function [v] = recrs_5(y0,y1,z)
r = zeros(1,z+1); %Vector to store all the values
r(1) = y0; %set first entry to y_0
r(2) = y1; %second entry to y_1
firstnegative = 0; %Finds first instance of negative value
for i = 3:z+1
    r(i) = ((5/2).*r(i-1) - r(i-2));
    if r(1,i)<=0 && firstnegative==0;
       i;  %spits out negative value occurance, * Note this is the reason why the recusrion goes to Infty.
       firstnegative=1;
    elseif r(i)>0;
        
    end
end
v = r;
end
