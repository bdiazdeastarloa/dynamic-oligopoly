% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
%
% main:
% - set parameters.
% - solve the model.
% - simulate and generate moments.
%
% Can be used standalone (i.e. not within estimation routine) to compute
% counterfactuals. See below for commented lines of code.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2014.
% =========================================================================

%% Settings and parameters.

% Calibrated parameter value (comment out when estimating).
parguess = [0.854664     0.472961     0.601777   0.00579473      3.39835      4.69322     0.987369      15.4188     0.869374];

% Assign new guess to parameter values.
alpha1  = parguess(1);              % r&d hazard: scale
alpha2  = parguess(2);              % r&d hazard: elasticity
eta1    = parguess(3);              % learning hazard: scale
eta2    = parguess(4);              % learning hazard: elasticity
etax1   = parguess(5);              % exit opportunity hazard
etax2   = 0;                        % exit opportunity hazard
etae1   = parguess(6);              % entry opportunity hazard
etae2   = 0 ;                       % entry opportunity hazard
delta   = parguess(7);              % negative productivity shock hazard
c0      = parguess(8);              % cost function level
beta    = parguess(9);              % cost function curvature
Phie_hi = 150;                      % entry cost upper bound
Phie_lo = 115;                      % entry cost lower bound
Phi_hi  = 100;                      % scrap value upper bound
Phi_lo  = 50;                       % scrap value lower bound
alpha   = 2.779;                    % price coefficient

% Set parameters, build state space (only if not estimating the model).
% Indicate scenario.
counter = 0;
% Assign parameters.
setpar;


%% Compute initial values for policies and value function.

initfile2 = ['initials_a' sprintf('%5.4f', par.alpha) '.mat'];
sol_last  = 'lastloop_sol.mat';

if exist(initfile2, 'file')==0
    % Set initial prices and values to static game solution
    [p0,V0] = init(par);
    x0 = zeros(par.N,par.S);
    y0 = zeros(par.N,par.S);
    eval(['save ',initfile2,' p0 V0 x0 y0']);
end

%if exist(sol_last,'file')
try
    % Load converged values from last parameter guess.
    load(sol_last)
    V0 = V1;
    p0 = p1;
    x0 = x1;
    y0 = y1;
    clear V1 p1 x1 y1
catch
%elseif exist(initfile2, 'file')
    % Load static equilibrium if it exists.
    load(initfile2)
end

% Assign initial values to par structure.
par.V0 = V0;
par.x0 = x0;
par.p0 = p0;
par.y0 = y0;

clear V0 x0 p0 y0

%% Solve for an equilibrium and simulate.

% Compile C file.
% MAC
% if exist('compMPE.mexmaci64','file')==0
%     mex compMPE.c -I/usr/local/Cellar/gsl/1.16/include -L/usr/local/Cellar/gsl/1.16/lib -lgsl -lgslcblas
% end
% LINUX (clusters)
% if exist('compMPE.mexa64','file')==0
%     mex file.c -I/usr/global/gsl/1.16/include -L/usr/global/gsl/1.16/lib -lgsl -lgslcblas
% end

% Solve for a MPE.
[info_,V1,x1,p1,y1] = compMPE(par);

solveflag = info_(1);
% Save converged results as initial condition for next guess.
if solveflag==1
    existence(par,V1,p1);
    save(sol_last,'V1','p1','x1','y1');
else
    disp('WARNING: solution algorithm did not converge.');
    disp('Punishment applied to moments.');
end

clear V1

% Simulate and compute moments.
[MM,~] = moms(par,p1,x1,y1,solveflag);

clear p1 x1 y1 par

% EOF
