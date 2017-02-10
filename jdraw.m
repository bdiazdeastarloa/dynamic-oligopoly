%% Get draw for jump (index of which event occurred).
function [draw,nj] = jdraw(x,prob)

% Replicates 'mnrnd(1,x)' to get random draw for new jump, but providing
% the random draw 'prob'.
% 'x' is a vector containing the probabilities of jumping to a specific
% event.

% x has to be a column vector. 
if size(x,2)~=1
    x = x';
end
% # of possible jumps.
nj   = size(x,1);                   
edge = [zeros(1,1);cumsum(x,1)];
draw = histc(prob,edge);
draw(:,end) = [];
draw = find(draw==1);
end