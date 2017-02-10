%
% Decode states (d,i1,i2,...,iN). 
%
function state = alldecode(code)

% Globals.
global D N M binom;

% Do it.
K = length(code);
state = zeros(N+1,K);
for k=1:K
    state(:,k) = decode(code(k));
end
