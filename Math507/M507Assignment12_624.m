% A is the matrix of T with respect to S given in the problem
A = sym([1 , 1, 1
    1 , 1 , 1
    1 , 1 , 1]);
% C will be an orthonormal basis for the range of A
C = orth(A)
% The output is C = 
% [3^(1/2)/3
% 3^(1/2)/3
% 3^(1/2)/3]
% which is a 3x1 column vector, indicating that the Rank(A) = 1
% This is not surprising because the eigenvalues of A are 3 and 0, where 0
% has a multiplicity of 2, therefore the nullity(A) = 2.
% Since the nullity(A) = 2, but T: R^3 to R^3, I do not believe that A is
% diagonalizable.