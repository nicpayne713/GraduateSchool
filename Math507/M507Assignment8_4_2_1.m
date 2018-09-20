% Problem 4.2.1a
% Matrix of the transformation
T = [1 -2 4
    -1 3 -6
    0 2 -3]
% vector z as defined in the problem
z = [1;0;0]
% finding 4 vectors to guarantee a linear dependence relation
z1 = T*z
z2 = T*z1
z3 = T*z2
% puts the vectors in a matrix to find the minimal polynomial
A = [z z1 z2 z3];
% prints the matrix in RREF to easily read off the coefficients of the
% minimal polynomial
Ar = rref(A)

%%
% Problem 4.2.1b
% same transformation matrix
T = [1 -2 4
    -1 3 -6
    0 2 -3]
% identity matrix
I = [1 0 0
    0 1 0
    0 0 1];
% vector u defined as in problem
u = (T*T+I)*z
% finding 3 other vectors along with u to guarantee a linear dependence
% relation
u1 = T*u
u2 = T*u1
u3 = T*u2

B = [u,u1,u2,u3];
% prints the matrix made up of column vectors resulting from applying T to
% u 3 times in RREF to easily read the coefficients for the minimal
% polynomial
Br = rref(B)