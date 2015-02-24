%function Y = main()

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
%
% Main file:
% - set parameters.
% - solve the model.
% - simulate and generate moments.
%
% Written by Bernardo Diaz de Astarloa @ PSU May-July 2014.
% =========================================================================

%% Settings and parameters.

clear all
clc

fprintf('Setting parameters...');
% Baseline scenario.
counter = 0;                    %#ok   

% Preliminary calibration values.
par0 = [0.01117   11.46020   0.79771   0.81111];
% Assign parameter values to solve for.
delta   = 0;               %#ok negative productivity shock
alpha1  = 0;               %#ok productivity hazard parameter
alpha2  = 0;               %#ok productivity hazard parameter
eta1    = par0(1);              %#ok productivity jump hazard
eta2    = 1;                    %#ok productivity jump hazard (0=no learning)
c0      = par0(2);              %#ok cost function level
etax1   = par0(3);              %#ok exit opportunity hazard
etae1   = par0(4);              %#ok entry opportunity hazard
Phi_hi  = 50;                   %#ok scrap value upper bound
Phi_lo  = 20;                   %#ok scrap value lower bound
Phie_hi = 100;                  %#ok entry cost upper bound
Phie_lo = 50;                   %#ok entry cost lower bound
alpha   = 2.779;                %#ok price coefficient

% Set parameters, build state space, set log files.
setpar;

%load('results/res_GS_N8M6D5Ly0.mat')
%V0 = V1;
%p0 = p1;

% Save variables and clear them from workspace.
filestr = sprintf('ind_N%dM%dD%d.mat',[N,M,D]);
eval(['save ',filestr,vars]);
eval(['clear ',vars]);

fprintf(' done.\n');

%% Solve for an equilibrium and simulate.

% Compile C file.
mex compMPE.c -I/usr/local/Cellar/gsl/1.16/include -L/usr/local/Cellar/gsl/1.16/lib -lgsl -lgslcblas

% Solve for a MPE.
compMPE(filestr);
load(filestr)

% Here I should check convergence in 'info' and assign 'solveflag'.
%{
% Compute moments.
MM = moms(par,P,X,E,solveflag);
if solveflag==0
    % Save converged results as initial condition for next guess
    save([codepath '/results/lastloop_sol.mat'],'P','V','X','E');
else
    disp('WARNING: there was a problem when solving the model.');
    disp('Punishment applied to moments.');
end

% Data and weights.
[data,W] = data_stats;

% Compute loss function
try
    error = data-MM;
    Y     = error'*W*error;
    Y_unw = norm(error)/norm(data);
    
    % Print Diagnostics
    compare = cat(2,data,model);
    
    fprintf('==============================================================\n');
    fprintf('Partial diagnostics:\n');
    fprintf('\n')
    fprintf('\nData - Model comparison:\n');
    disp(   num2str(compare));   
    fprintf('\nParameters values:\n');
    disp(   num2str(par0));
    fprintf('\nLoss (weighted):\n');
    disp(   num2str(Y));
    fprintf('\nLoss (unweighted):\n');
    disp(   num2str(Y_unw));
    fprintf('==============================================================\n');

catch ER 

    % Something went wrong: high loss
    disp('Something went wrong when computing metric: high loss.');
    Y = 1e10;    
end
 %}
