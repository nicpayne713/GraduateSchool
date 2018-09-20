%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART a)
f = @(x) exp(x)/x;
truea = expint(-1)-expint(-2);
ra = romberg(f, 1, 2, 5);
Ea = ra-truea;
% This formats the error table so that the upper entries are all 0
for i = 1:5;
    for j = 1:5;
        if i<j;
            Ea(i,j)=0;
        end
    end
end
ra % prints the romberg table and error table for part a
Ea
%
% PART b)
g = @(x) sqrt(1-x.^2);
trueb = integral(g, 0,1);
rb = romberg(g, 0,1,5);
Eb = rb-trueb;
for i = 1:5;
    for j = 1:5;
        if i<j;
            Eb(i,j)=0;
        end
    end
end
rb % prints the romberg table and error table for part b
Eb