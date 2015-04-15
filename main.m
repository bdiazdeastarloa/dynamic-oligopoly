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
parguess = [0.131   1.0900   2.500   2.500   11.1   0.7];

% Assign new guess to parameter values.
alpha2  = parguess(1);              % r&d hazard: elasticity
alpha1  = 1;                        % r&d hazard: scale
eta1    = parguess(2);              % learning hazard: scale
eta2    = 0;                        % learning hazard: elasticity
etax1   = parguess(3);              % exit opportunity hazard
etax2   = 0;                        % exit opportunity hazard
etae1   = parguess(4);              % entry opportunity hazard
etae2   = 0 ;                       % entry opportunity hazard
c0      = parguess(5);              % cost function level
delta   = parguess(6);              % negative productivity shock hazard
Phi_hi  = 100;                      % scrap value upper bound
Phi_lo  = 50;                       % scrap value lower bound
Phie_hi = 150;                      % entry cost upper bound
Phie_lo = 100;                      % entry cost lower bound
alpha   = 2.779;                    % price coefficient

% Set parameters, build state space (only if not estimating the model).
% Indicate scenario.
counter = 0;
% Assign parameters.
setpar;

% Compute initial values for policies and value function.
initfile2 = ['initials_N' int2str(N) 'M' int2str(M) 'D' int2str(D) 'a' sprintf('%5.4f', alpha) '.mat'];
sol_last  = ['lastloop_N' int2str(N) 'M' int2str(M) 'D' int2str(D) '.mat'];

if exist(sol_last,'file')
    % Load converged values from last parameter guess.
    load(sol_last)
    V0 = V1;
    p0 = p1;
    x0 = x1;
    y0 = y1;
    clear V1 p1 x1 y1 info_
elseif exist(initfile2, 'file')
    % Load static equilibrium if it exists.
    load(initfile2)
else
    % Set initial prices and values to static game solution
    [p0,V0] = init(par);
    x0 = zeros(N,S);
    y0 = zeros(N,S);
    eval(['save ',initfile2,' p0 V0 x0 y0']);
end

% Assign initial values to par structure.
par.V0 = V0;
par.x0 = x0;
par.p0 = p0;
par.y0 = y0;

% Save variables and clear them from workspace.
% eval(['save ',mexinp,vars]);
% eval(['clear ',vars]);


%% Solve for an equilibrium and simulate.

% Compile C file.
if exist('compMPE.mexmaci64','file')==0
    % mex compMPE.c -I/usr/local/Cellar/gsl/1.16/include -L/usr/local/Cellar/gsl/1.16/lib -lgsl -lgslcblas
    mex compMPE.c -I/usr/local/include -L/usr/local/lib -lnlopt -lm
end

% Solve for a MPE.
[info_,V1,x1,p1,y1] = compMPE(par);

solveflag = info_(1);
% Save converged results as initial condition for next guess.
if solveflag==1
    save(sol_last,'V1','p1','x1','y1','info_');
else
    disp('WARNING: solution algorithm did not converge.');
    disp('Punishment applied to moments.');
end

% Simulate and compute moments.
MM = moms(par,p1,x1,y1,solveflag);

% EOF