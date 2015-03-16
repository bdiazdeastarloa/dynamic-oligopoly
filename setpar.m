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
ie0 = 3;                                         % productivity level at entry (fixed case)

% Estimate import price process.
% gamma: hazard rate of a jump. 
% p0vec: grid.
% Q: intensity matrix.
load p0;                                         % load imported price series
[~,p0vec,Qp] = p0process(p0,2);                  % estimate imported price process 
% p0vec  = 2.3;
% Qp     = 0;
% load p0process_alt;
% p0vec = p0vec_alt;
% Qp    = Qp_alt;
p0size = size(p0vec,1);                          % # of imported price states 
if counter==1                                    % countefactual 1:
    tariff = 0.3;                                % 30% tariffs
    p0vec  = (1+tariff)*p0vec;     
end
clear p0;

% Demand.
% mkt = [100 160 220];                          
mkt = 160;
dsize = size(mkt,2);                             % # of demand states     
Qd = zeros(dsize,dsize);                         % demand transitions
if dsize>1
    Qd = zeros(dsize,dsize)-log(2/3).*diag(ones(dsize-1,1),1);
    Qd(dsize,1) = -log(2/3);
end
D = dsize*p0size;                                % # of (d,p0) common states

% Common states transition matrix (DxD).
trans = kron(Qd,eye(p0size,p0size));
for i=1:dsize
    i1 = 1+(i-1)*p0size;
    i2 = i1+p0size-1;
    trans(i1:i2,i1:i2) = Qp;
end
% Hazard of an aggregate shock (sum of a row of 'trans').
gamma = sum(trans,2);
gamma = gamma(1);

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
lambday = 0.9;                                   % dampening factor: exit/entry policies   
maxiter = 500;                                   % max # of iterations
tol     = 1e-4;                                  % tolerance
steps = 0;                                       % # steps in modified policy iteration.


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
% (This is the cartesian product of the indices.)
[agg1,agg2] = meshgrid((1:dsize),(1:p0size));
agg = [agg1(:) agg2(:)];


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
        ' etax2 ',...                   % Exit: hazard rate param.    
        ' etae1 ', ...				    % Entry: hazard rate param.
        ' etae2 ', ...                  % Entry: hazard rate param.
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

% File name to give as input to MEX-file.
mexinp = sprintf('ind_N%dM%dD%d.mat',[N,M,D]);


%% Create parameter structure to pass to simulations.

par = struct();
par.w      = w;
par.rho    = rho;
par.alpha  = alpha;
par.N      = N;
par.M      = M;
par.D      = D;
par.S      = S;
par.p0size = p0size;
par.dsize  = dsize;
par.mkt    = mkt;
par.p0vec  = p0vec;
par.mgc    = mgc;
par.trans  = trans;
par.binom  = double(binom);
par.state  = double(state);
par.agg    = agg;
par.delta  = delta;
par.alpha1 = alpha1;
par.alpha2 = alpha2;
par.eta1   = eta1;
par.eta2   = eta2;
par.etax1  = etax1;
par.etax2  = etax2;
par.etae1  = etae1;
par.etae2  = etae2;
par.gamma  = gamma;
par.ie0    = ie0;
par.Phi_hi = Phi_hi;
par.Phi_lo = Phi_lo;
par.Phie_hi= Phie_hi;
par.Phie_lo= Phie_lo;