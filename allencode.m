%
% Encode states (d,i1,i2,...,iN). 
%
function [code,pos] = allencode(state,binom,N,D)

% Do it.
K = size(state,2);
code = zeros(1,K);
pos = zeros(N,K);
for k=1:K
    [code(k),pos(:,k)] = encode(state(:,k),binom,N,D);
end
