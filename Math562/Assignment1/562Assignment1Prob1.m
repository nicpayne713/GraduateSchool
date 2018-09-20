% Math 562 Assignment 1 Problem 1
a = [2  0 0 0 ; 0 1 0 0 ; 0  0  1 0 ;  0  0 0 1   ];% double column 1
b = [1  0 0 0 ; 0 1 0 0 ; 0  0 .5 0 ;  0  0 0 1   ];% half row 3
c = [1  0 1 0 ; 0 1 0 0 ; 0  0  1 0 ;  0  0 0 1   ];% add row 3 to row 1
d = [0  0 0 1 ; 0 1 0 0 ; 0  0  1 0 ;  1  0 0 0   ];% switch columns 1 & 4
e = [1 -1 0 0 ; 0 1 0 0 ; 0 -1  1 0 ;  0 -1 0 1   ];% - row 2 from others
f = [1  0 0 0 ; 0 1 0 0 ; 0  0  1 1 ;  0  0 0 0   ];% replace column 4 by 3
g = [0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ];               % delete column 1

B = [1 2 3 4;5 6 7 8;2 4 6 8;1 3 5 7]; % Some random B to check answers

e*c*b*B*a*d*f*g; % order the matrices appear in with row operations on the 
                 % left of B and column operations on the right
A = e*c*b;
C = a*d*f*g;