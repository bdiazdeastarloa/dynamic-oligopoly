%% Compute consumer surplus.
function cs = surplus(par,p,p0,d)

cs = -d*(1/par.alpha)*log(sum(exp(-par.alpha*p))+exp(-par.alpha*p0));

end