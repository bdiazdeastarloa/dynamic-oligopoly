%
% Decode state (d,i1,i2,...,iN).
%
function state = decode(code,binom,N,D)

% Compute subscript (d,j) in a Matlab array of size (D,J) from linear index.
state = ones(N+1,1);
code = code-1;
j = floor(code./D);
state(1) = rem(code,D)+1;

% Map j=0,...,J-1 into (i1,i2,...,iN).
for n=N+1:-1:2
    while (binom(state(n)+n-1,state(n)+1)<=j)
        state(n) = state(n)+1;
    end
    j = j-binom(state(n)+n-2,state(n));
end
