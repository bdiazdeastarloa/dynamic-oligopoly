% =========================================================================
% Dynamic oligopoly with learning by doing
% Continuous time version
%
% Calibration routine (genetic algorithm)
% 
% - calls 'distance_lbd', which contains main solution routine and moment
%   generating routines.
%
% Written by Bernardo Diaz de Astarloa @ PSU September 2014.
% =========================================================================

clear all;
close all;
clc;

tic;
filename = 'results_ga-estimation.mat';
format long;

rng(82082);

fprintf('==============================================================\n')
fprintf('Dynamic oligopoly industry dynamics.\n');
fprintf('Continuos time.\n');
fprintf('\nGenetic algorithm estimation.\n');
fprintf('==============================================================\n')

% Initial condition (# of individuals x size of parameter vector).
% [rd_par1 rd_par2 lbd_par1 lbd_par2 etax etae delta c0 beta].

pop = [];

nvar = size(pop,2);
% Lower and upper bounds for search.
lbound = [0.5000   0.30000   0.4000   0.00200   2.00000   2.0000    0.4   10.000    0.6]; 
ubound = [1.0000   0.6000   0.90000   0.00600   8.0000   8.000   1   16.000   1.0];                 

% Set GA algorithm options
% Stop if no improvement for (secs)...
stall = 3*3600;
% Kill algorithm after (secs)...
kill = 24*3600;

options = gaoptimset('Display','iter','PopulationSize',8,'Generations',300,...
'StallTimeLimit',stall,'TimeLimit',kill,'MutationFcn',@mutationadaptfeasible,...
'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
'PlotFcns',@gaplotbestf,'EliteCount',0,'HybridFcn',@fmincon);

[Y,fval,exitflag,output,population,scores] = ga(@(Y) distance(Y),nvar,[],[],[],[],lbound,ubound,[],options);  

ga_time = toc/3600;
disp(['The job took ',num2str(ga_time),' hours.']);
eval(['save ',filename,' Y fval exitflag output population scores ga_time']);
