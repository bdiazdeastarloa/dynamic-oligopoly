function f = foc_static(p,pars)

% =========================================================================
% First order conditions for prices
%
% Used to solve for optimal prices, given an industry structure and an 
% expected continuation value.
% 
% Inputs:
% - p is the initial price vector.
% - par is the parameter structure.
% - active is the index of active firms to solve for.
% - w is the industry structure.
% - deltav is the increase in the value upon an increase in the state.
% 
% Written by Bernardo Diaz de Astarloa @ PSU June 2014
% =========================================================================

a  = pars.a;
p0 = pars.p0;
d  = exp(-a*(p-p0));               
s  = d./(1+sum(d));
c  = pars.c;

f = 1 - a.*(1-s).*(p-c);