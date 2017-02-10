function f = foc_j(pj,j,pp,deltav)

% =========================================================================
% First order conditions for prices
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

p     = pp.p;
p(j)  = pj;
alpha = pp.alpha;                 
nu    = pp.nu;
eta   = pp.eta;
p0    = pp.impp;
denom = 1+sum(exp(- alpha*(p-p0)));               
s     = exp(- alpha*(p-p0))./denom;
sj    = s(j);
mgc   = pp.c(j);
q     = pp.dem*s;
deltavj = deltav(j);              % value of increase in own state 

q1    = eta*nu*q.^(nu-1);         % derivative of hazard function
q1j   = q1(j);
q2    = q1.*s;

f = (1/alpha) - (1-sj)*(pj-mgc) - q1j*deltavj + q2'*deltav;

