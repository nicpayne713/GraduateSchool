% Problem 2.6.2 on page 84
% A is the matrix of the transformation S(a+bx+cx^2) = (a+2b+c) +
% (2a+3b+2c)x + (a+3b+2c)x^2
%
% Take B to be the standard basis for R_2 [x] = {1, x, x^2}
%S(1) = 1 + 2x + x^2
%S(x) = 2 + 3x + 2x^2
%S(x^2) = 1 + 3x + 2x^2
%
% So then A is:
A = [1,2,1
2,3,2
1,3,2];

% if S is invertible, then A will be also
B = inv(A);

% So B is the inverse of A, therefore the transformation that B describes
% will be S^-1
%
% To do this multiply by the vector [a b c]^T
syms a b c;
v = [a,b,c]';
C = B*v

% After running the code we have that the matrix C is
% C = [0  1 -1
%      2 -1  0
%     -3  1  1]


