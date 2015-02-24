% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
%
% Settings and parameters.
%
% Written by Bernardo Diaz de Astarloa @ PSU May 2014.
% Adapted from code by Gowrisankaran (1993) edited by Ron Goettler (2011).
% ========================================================================


%% Main parameters.

N   = 8;                                         % max # of firms  
rho = -log(0.925);                               % discount rate

% State variables
w = [1 2.3 3.6 5 8.1 11.5];                      % productivity vector
M = size(w,2);                                   % # of productivity states
% M  = 6;
% wl = 1;
% wh = 11.5;
% w  = [wl:(wh-wl)./(M-1):wh]';
ie0 = M-2;                                       % productivity level at entry (fixed case)

% Estimate import price process.
% gamma: hazard rate of a jump. 
% p0vec: grid.
% Q: intensity matrix.
load p0;                                         % load imported price series
[gamma,p0vec,Qp] = p0process(p0,2);              % estimate imported price process 
if counter==1                                    % countefactual 1:
    tariff = 0.3;                                % 30% tariffs
    p0vec  = (1+tariff)*p0vec;     
end
% p0vec  = 2.3;
% Qp     = 0;
p0size = size(p0vec,1);                          % # of imported price states 

% Demand.
% mkt = [100 160 220];                          % domestic demand (MW)
mkt = 160;
dsize = size(mkt,2);                             % # of demand states     
Qd = zeros(dsize,dsize);                         % demand transitions
if dsize>1
    Qd = zeros(dsize,dsize)-log(2/3).*diag(ones(dsize-1,1),1);
    Qd(dsize,1) = -log(2/3);
end
D = dsize*p0size;                               % # of (d,p0) common states

% Common states transition matrix (DxD).
trans = kron(Qd,eye(p0size,p0size));
for i=1:dsize
    i1 = 1+(i-1)*p0size;
    i2 = i1+p0size-1;
    trans(i1:i2,i1:i2) = Qp;
end

% Technology.
beta = -1;                                       % slope of learning curve
if counter==2                                    % counterfactual 2:                    
    c0 = 0.8*c0;                                 % 20% reduction in marginal cost
end
mgc = (c0*w.^beta)';                             % marginal cost schedule

% Contraction routine settings
method = 'GaussSeidel';
lambdaV = 0.9;                                   % dampening factor: value function
lambdax = 0.9;                                   % dampening factor: investment policy
lambdap = 0.9;                                   % dampening factor: pricing policy   
lambday = 0.8;                                     % dampening factor: exit/entry policies   
maxiter = 200;                                   % max # of iterations
tol     = 1e-4;                                  % tolerance
steps = 0;                                       % # steps in modified policy iteration.

% Simulation settings
T    = 11;                                       % number of periods to simulate
sims = 500;                                      % number of simulations

% Random draws for simulations (T*30,sims)
% if exist('rand_draws.mat','file')
%     load rand_draws;
% else
%     rng(8282);
%     jumps     = 30;                              % allow for many jumps per period
%     draws_agg = rand(T*jumps,sims);              % draws: aggregate jump 
%     draws_ind = rand(T*jumps,sims);              % draws: firm-level jumps
%     draws_out = mnrnd(1,[0.5 0 0.5],330);        % draws: outside good     
%     save rand_draws draws_agg draws_ind draws_out;
% end
d0  = 1;                                         % initial demand state
p00 = 1;                                         % initial foreign price state
w0  = [1 1 1 2 4 M+1 M+1 M+1];                   % initial industry structure
common0 = [d0 p00];                              % initial common state


%% State space encoding.

filename = ['encoding_N' int2str(N) 'M' int2str(M) 'D' int2str(D) '.mat'];
if exist(filename,'file')
    load(filename)
else
    % Total # of unique industry structures = D*[(M+1+N-1)!/(M+1-1)!N!] (P&M 1994).
    S = D*( factorial(M+1+N-1)/(factorial(M+1-1)*factorial(N)) );

    % Setup binomial coefficients for decoding/encoding of states.
    binom = [zeros(N+M+1,1),eye(N+M+1)];
    n = 2;
    while (n<=N+M+1)
        binom(n,2:n) = binom(n-1,2:n)+binom(n-1,1:n-1);
        n = n+1;
    end

    % Decode state indices into industry structures table.
    state = zeros(N+1,S);
    for s=1:S
        state(:,s) = decode(s,binom,N,D);
    end
        
    % Convert to unsigned 64-bit integer (ulong).
    % NOTE: for Windows the compatible type is 32-bit.
    binom = uint64(binom);
    % Convert to unsigned 8-bit integer (uchar).
    state = uint8(state);
    
    eval(['save ',filename,' S binom state']);
end

% Demand/Foreign price combinations indices.
agg = zeros(D,2);
for i=1:dsize
    for j=1:p0size
        ind = j + (i-1)*p0size;
        agg(ind,1) = i;
        agg(ind,2) = j;
    end
end

%% Compute initial values for policies and value function.

initfile2 = ['initials_N' int2str(N) 'M' int2str(M) 'D' int2str(D) 'a' sprintf('%5.4f', alpha) '.mat'];
if exist(initfile2, 'file')
    load(initfile2);
else
    state_ = double(state);
    % Set initial prices and values to static game solution
    [p0,V0] = init(N,M,S,state_,agg,mgc,p0vec,mkt,alpha,rho);
    eval(['save ',initfile2,' p0 V0']);
    clear state_
end

x0 = zeros(N,S);
y0 = zeros(N,S);

%% Declare variables names.

% Declare variable names to pass to MEX-file.
% NOTE: This is a vector of strings, not mat files.

vars = [
        ' rho ', ...					% Discount factor.
        ' alpha ',...                   % Demand system: price coefficient.
        ' N ',...                       % # of firms.
        ' M ',...                       % # of productivity states.
        ' D ',...                       % # of common states.
        ' S ',...                       % # of unique industry states.
        ' p0size ', ...					% Common state: # of foreign price values.
        ' dsize ', ...					% Common state: # of domestic demand values.
        ' mkt ',...                     % Common state: domestic demand values.
        ' p0vec ',...                   % Common state: foreign price values.
        ' mgc ',...                     % Marginal cost values (for each productivity level).
        ' trans ', ...					% Transition rates: common state.
        ' binom ',...                   % State space encoding: binomial table.
        ' state ',...                   % State space encoding: list of ind. structures.
        ' delta ', ...					% Productivity process: negative shock.
        ' alpha1 ', ...					% Productivity process: hazard rate param.
        ' alpha2 ', ...					% Productivity process: hazard rate param.
        ' eta1 ', ...					% LBD process: hazard rate param.
        ' eta2 ', ...					% LBD process: hazard rate param.
        ' etax1 ', ...					% Exit: hazard rate param.
        ' etae1 ', ...				    % Entry: hazard rate param.
        ' ie0 ', ...					% Potential entrant initial state.
        ' Phi_hi ', ...					% Scrap value: upper bound.
        ' Phi_lo ', ...					% Scrap value: lower bound.
        ' Phie_hi ', ...				% Setup cost: upper bound.
        ' Phie_lo ',...                 % Setup cost: lower bound.
        ' method ',...                  % Program control: Method.
        ' maxiter ', ...                % Program control: Maximum number of iterations.
        ' tol ', ...                    % Program control: Tolerance.
        ' steps ', ...				    % Program control: # steps in modified policy iteration.
        ' lambdaV ', ...                % Program control: Weight of updated value function.
        ' lambdax ', ...                % Program control: Weight of updated investment policy.
        ' lambdap ', ...                % Program control: Weight of updated pricing policy.
        ' lambday ', ...                % Program control: Weight of updated entry/exit policy.
        ' V0 ', ...						% Starting values: Value function.
        ' x0 ', ...						% Starting values: Investment policy.
        ' p0 ',...                      % Starting values: Pricing policy.
        ' y0 '  					    % Starting values: Entry/exit policy.
        ];
    
 