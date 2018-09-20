% 561 Final number 2
%b
syms h lambda
I = [1,0;0,1];
c = [1/2;1];
c2 = [0;1];
DG = [h*lambda,-1];
B = [1,1;0,1];

%PEC
Spec = (I + c*DG)*B

%PECE
Spece = (I + c2*DG)*(I + c*DG)*B

% check for zero stability by using h*lambda = 0
subs(Spec,h*lambda,0)
subs(Spece,h*lambda,0)

% Spec has eigen values of 1 and -1, whose magnitidues are less or equal 1
% therefore it is zero stable, and Spece has eigen values 1 and 0, whose
% magnitudes are less than 1 therefore it's also zero stable.

% To see how far along the negative real axis the regions of stability
% stretch we use the following code
syms x % symbolic representation of hlambda used to find eigenvalues of the 
% S matrices which govern the regions of stability
dg = [x ,-1]; %symbolic DG matrix
S1 = (I + c*dg)*B;
S2 = (I + c2*dg)*(I + c*dg)*B;
lamb1 = 0;
lamb2 = 0;
z=[];
while isempty(find(abs(double(z))>1)) % only can have roots less or equal 1
z = eig(subs(S1,lamb1));
lamb1 = lamb1 - .005; % accounts for root which was barely greater than 1
end
z=[];
while isempty(find(abs(double(z))>1))
z = eig(subs(S2,lamb2));
lamb2 = lamb2 - .005;
end
lamb1,lamb2