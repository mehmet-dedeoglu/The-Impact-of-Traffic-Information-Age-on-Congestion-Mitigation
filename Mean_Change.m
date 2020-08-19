function change = Mean_Change(age,dT,DAR,ell)



flatNormAge = floor(age/dT);
change = (age/dT-flatNormAge)*DAR(ell-flatNormAge)+...
    sum(DAR(ell-flatNormAge+1:ell+1));

end