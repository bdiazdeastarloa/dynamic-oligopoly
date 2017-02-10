% =========================================================================
% Dynamic oligopoly with learning by doing
% Continuous time version
%
% Compute standard errors using bootstrapped covariance matrix of moments.
%
% Written by Bernardo Diaz de Astarloa @ PSU April 2015.
% =========================================================================

clc
tic
% Load bootstrapped covariance matrix and associated parameters.
% load results_ga-estimation
% clearvars -EXCEPT Y       % Keep parameter vector only.
Y = [0.727799     0.466028     0.586318   0.00337886      4.08584      4.08302     0.598165      12.1315     0.783536];
load bootcov

% Parameters. 
S = 500;
nummoms = size(bootcov,1);
pop = Y; 
findif = 1e-3;
jac = zeros(nummoms,size(pop,2));

% Weighting matrix and moments difference.
[~,W,base] = distance(pop);

% Compute derivatives of moments difference.
basis = eye(size(pop,2));
for k = 1:size(pop,2)
    display(k);
    ptemp = pop .* (1 + basis(k,:) * findif);
    [~,~,new] = distance(ptemp);
    jac(:,k) = (new-base) / (pop(k) * findif);
end

% Compute covariance matrix.
varcov = (1+1/S)*(jac' * W * jac)^-1 * (jac' * W * bootcov * W * jac) * (jac' * W * jac)^-1;
% Use this if using the optimal weighting matrix W = bootcov^(-1):
% varcov = (1+1/S)*(jac' * W * jac)^-1;
posdefchk = min(eig((varcov+varcov')/2));
if posdefchk<0
    display('WARNING: Not a local maximum.');
end
se = sqrt(diag(varcov));

toc
save 'results/std_errors_20161201_tight.mat'
