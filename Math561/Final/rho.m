function r = rho(theta,r0)

% For characteristic polynomial CPOLY (which is a function of Z)
% set Z = exp(i*THETA), and find the value R for which the largest root
% of cpoly exceeds 1 for the first time. R0 is an initial guess.
%
% This gives one point in the polar coordinate form for the boundary
% of the region of stability.

global cpoly z
if (nargin < 2)
    r0 = 1;
end

c = @(t) max(abs(roots(double(subs(cpoly,z,exp(1i*theta)*t)))))-1;
% The following IF statement keeps FZERO from wandering off
% towards infinity in the case THETA = pi/2.
if (theta == pi/2)
    try
	r = fzero(c,[r0/2,2*r0]);
    catch
	r = 0;
    end
else
    r = fzero(c,r0);
end