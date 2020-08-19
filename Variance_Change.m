function [change,changeEst] = Variance_Change(age,dT,DVAR,ell)

flatNormAge = floor(age/dT);
change = (age/dT-flatNormAge)*DVAR(ell-flatNormAge)+...
    sum(DVAR(ell-flatNormAge+1:ell+1));

changeEst = (age/dT-flatNormAge)*DVAR(ell-flatNormAge)+...
    sum(DVAR(ell-flatNormAge+1:ell));

end