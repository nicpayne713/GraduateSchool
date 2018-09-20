function f = alpha1D(uplus, uminus)
% function for alpha parameter in 1D LLF numericl Hamiltonian from
% Osher_Shu_ENO_Schemes
% Input is the local points 
% For now, this needs to be updated with an exact derivative fo the
% Hamiltonian being used.

    dH = @(u) u+1;

    x = min(uplus, uminus);
    y = max(uplus, uminus);
    I = linspace(x,y,1000);
    Hdiff = abs(dH(I));
    f = max(Hdiff);