%
% Encode state (d,i1,i2,...,iN). 
%
function [code,pos] = encode(state,binom,N,D)

% Map (i1,i2,...,iN) into j=0,...,J-1.
j = 0;
[state(2:N+1),pos] = sort(state(2:N+1)); 
pos(pos) = (1:N)';
for n=N+1:-1:2
    j = j+binom(state(n)+n-2,state(n));
end

% Compute linear index from subscript (d,j) in a Matlab array of size (D,J).
code = 1+(state(1)-1)+D.*j;
