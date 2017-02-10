function [gammahaz,p0vec,Q,jump] = pforprocess_ehr(data,n)

% =========================================================================
% p0process: estimates a discretized Ehrenfest diffusion process from
% discrete 'data' and a 2*n+1 grid size.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2014.
% =========================================================================

%if exist('p0process.mat', 'file')
%    load p0process;
%else
    [gammahaz,p0vec,Q,jump] = poisson_par(data,n);
    save p0process gammahaz Q p0vec
%end

end

function [gamma,p0vec,Q,jumpsize] = poisson_par(data,n)

delta = 1;
[mu,sig,lambda] = ou_calib(data,delta);

gamma    = lambda*n;        % arrival rate of shock
jumpsize = sig*gamma^(-.5); % size of jump states

% Build intensity matrix
Q = zeros(2*n+1);
Q(1,1) = -gamma;
Q(1,2) = gamma;
Q(2*n+1,2*n) = gamma;
Q(2*n+1,2*n+1) = -gamma;
for k=2:2*n
   Q(k,k) = -gamma;
   Q(k,k-1) = gamma*(.5*(1+((k-n-1)/(n))));
   Q(k,k+1) = gamma*(.5*(1-((k-n-1)/(n))));
end

% Zeros in the diagonal
temp = (Q<0);
Q(temp) = 0;

p0vec = mu + ((1:2*n+1)'-(n+1))*jumpsize;

end


function [mu,sig,lambda] = ou_calib(data,delta)
  n = length(data)-1;
 
  Sx  = sum( data(1:end-1) );
  Sy  = sum( data(2:end) );
  Sxx = sum( data(1:end-1).^2 );
  Sxy = sum( data(1:end-1).*data(2:end) );
  Syy = sum( data(2:end).^2 );
 
  mu  = (Sy*Sxx - Sx*Sxy) / ( n*(Sxx - Sxy) - (Sx^2 - Sx*Sy) );
  lambda = -log( (Sxy - mu*Sx - mu*Sy + n*mu^2) / (Sxx -2*mu*Sx + n*mu^2) ) / delta;
  a = exp(-lambda*delta);
  sigh2 = (Syy - 2*a*Sxy + a^2*Sxx - 2*mu*(1-a)*(Sy - a*Sx) + n*mu^2*(1-a)^2)/n;
  sig = sqrt(sigh2*2*lambda/(1-a^2));
  
end