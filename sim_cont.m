function out = sim_cont(par,P,X,Y,di,drw1,drw2,drw3,T)

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
% 
% sim_cont:
% - simulates arrival times and state jumps of the continuous time model.
% - inputs are: parameter structure, # of periods (T), policy 
%   functions (P, X, Y), and the initial state (w0,common0). 
% 
% Written by Bernardo Diaz de Astarloa @ PSU 2014
% =========================================================================

N = par.N;
M = par.M;
D = par.D;
binom = double(par.binom);
p0size = par.p0size;
p0vec  = par.pfor;
mkt    = par.mkt;
mgc    = par.mgc;
Q      = par.trans/par.gamma;      % aggregate state transition matrix
F = N*2;                           % max number of firms in history
TT = T*20;                         % max number of jumps in history

times_c = zeros(TT,1);             % arrival times
fprod_c = zeros(TT,F);             % firms' productivities
p_c     = zeros(TT,F);             % firms' prices
x_c     = zeros(TT,F);             % firms' investment
exit_c  = zeros(TT,1);             % exit events
entry_c = zeros(TT,1);             % entry events
p0_c    = zeros(TT,1);             % outside good price
q0_c    = zeros(TT,1);             % outside good quantity
d_c     = zeros(TT,1);             % outside good price
q_c     = zeros(TT,F);             % firms' quantities
rev_c   = zeros(TT,F);             % f irms' revenues
mgc_c   = zeros(TT,F);             % firms' marginal cost
pij_c   = zeros(TT,F);             % firms' profits
ecostj_c = zeros(TT,F);        
scrapj_c = zeros(TT,F);
ecost_c  = zeros(TT,1);            % entry costs
scrap_c  = zeros(TT,1);            % exit scrap values
cs_c     = zeros(TT,1);            % consumer surplus
ps_c     = zeros(TT,1);            % aggregate profits
nf_c     = zeros(TT,1);            % # active firms
rdshock  = zeros(TT,F);            % keep track of r&d shocks
lbdshock = zeros(TT,F);            % keep track of lbd shocks
tentry   = NaN(TT,F);            % record entry times

% Set initial state.
in  = (di(2:end)<M+1);             % index of incumbents
j   = sum(in);                     % # active firms
w   = encode(di,binom,N,D);       
d0  = di(1);
d   = ceil(d0/p0size);
p0  = d0 - (d-1)*p0size;
d   = mkt(d);
p0  = p0vec(p0);
t   = 0;

% Assign policies.
p   = P(:,w);
x   = X(:,w);
y   = Y(:,w);

% Compute profits, quantities and surplus.
[pr,q,q0] = profits(par,di(2:sum(in)+1),p(in),p0,d);
cs        = surplus(par,p(in),p0,d);
rev       = p(in).*q;
ps        = sum(pr);

% Compute hazards of jumps.
[hj,sh] = hazz(par,di,q,x,y);

% Initialize jumps #.
n = 1;                            
nextdi = di;

% Firms' positions.
id = zeros(N,1);
id(1:j) = (1:1:j);
% Start active firms counter (all firms ever present).
nf = j;

% Record initial state.
cs_c(n) = cs;                
ps_c(n) = ps;                
nf_c(n) = j;                
p0_c(n) = p0;
q0_c(n) = q0;
d_c(n)  = d;
nf_c(n) = j;

% Store firm-specific values.
fprod_c(n,id(1:j)) = nextdi(2:j+1);                     % productivity
p_c(n,id(1:j)) = p(1:j);                                % price
x_c(n,id(1:j)) = x(1:j);                                % rd investment
q_c(n,id(1:j)) = q(1:j);                                % quantitiy sold
rev_c(n,id(1:j)) = rev(1:j);
pij_c(n,id(1:j)) = pr(1:j);                             % profits
mgc_c(n,id(1:j)) = mgc(nextdi(2:j+1));                  % marginal cost
tentry(:,id(1:j))= t;      

flag = 0; 
n = n + 1;
t = - log(drw1(n))/sh;                           % first actual jump
if t==Inf
    disp('WARNING: t=Inf - no first jump. Something is wrong.');
    flag = 1;
end
   
while t>0 && t<T
    % There's a jump!
    times_c(n) = t;                   
    % Where to jump?
    % Get transition probabilities conditional on a jump.
    hh = hj./sh;
    [draw,nj] = jdraw(hh,drw2(n));
    % Search for specific shock going backwards.
    if draw==nj
        % Aggregate shock (demand or foreign price).
        trans = Q(d0,:);
        [d0,~] = jdraw(trans,drw3(n));
        nextdi(1) = d0;
    elseif draw==nj-1
        % Entry event.             
        we = max(nextdi(j+1)-2,1);         % entry draw
        %we = par.ie0; 
        nextdi(j+2) = we;
        nf = nf + 1;                       % add one firm to total firms
        id(j+1) = nf;                      % new firm ID
        % Update entry event and entry cost.
        entry_c(n) = 1;
        ecost_c(n) = par.Phie_lo + drw2(n)*(par.Phie_hi-par.Phie_lo);
        ecostj_c(n,id(j+1)) = ecost_c(n);
        tentry(n:end,id(j+1)) = t;
    elseif draw>3*j                        % 1st 2j cells are productivity jumps
        % Exit event.
        ii = draw-3*j;                     % which firm exits
        nextdi(ii+1) = M+1; 
        % Update exit event and exit cost.
        exit_c(n) = 1;
        scrap_c(n) = par.Phi_lo + drw2(n)*(par.Phi_hi-par.Phi_lo);
        scrapj_c(n,id(ii)) = scrap_c(n);
        id(ii) = 0;                         % free ID slot
    elseif draw>2*j
        % Negative productivity shock.
        ii = draw-2*j;                      % which firm goes down
        nextdi(ii+1) = nextdi(ii+1)-1;  
    elseif draw>j
        % R&D shock.
        ii = draw-j;                        % which firm goes up
        nextdi(ii+1) = nextdi(ii+1)+1; 
        % Update r&d shock tracker.
        rdshock(n,id(ii)) = 1;
    else %draw<=j
        % LBD shock.
        ii = draw;                          % which firm goes up
        nextdi(ii+1) = nextdi(ii+1)+1;
        % Update lbd shock tracker.
        lbdshock(n,id(ii)) = 1;        
    end
    
    % Sort new state.
    temp   = sortrows([nextdi(2:end) id]);    
    nextdi = [nextdi(1);temp(:,1)];
    id     = temp(:,2);              
       
    % Update everything.
    in = (nextdi(2:end)<M+1);  
    j  = sum(in);
    if j==0
        disp('WARNING: empty industry.');
        flag=1;
        break
    end
    w  = encode(nextdi,binom,N,D);
    p  = P(:,w);
    x  = X(:,w);
    y  = Y(:,w);
    d  = ceil(d0/p0size);
    p0 = d0 - (d-1)*p0size;
    d  = mkt(d);
    p0 = p0vec(p0);

    [pr,q,q0] = profits(par,nextdi(2:j+1),p(in),p0,d);
    cs        = surplus(par,p(in),p0,d);
    rev       = p(in).*q;
    ps        = sum(pr);

    [hj,sh] = hazz(par,nextdi,q,x,y);
    
    % Record aggregate outcomes.
    cs_c(n) = cs;                      % consumer surplus
    ps_c(n) = ps;                      % aggregate profits
    nf_c(n) = j;                       % # of active firms
    p0_c(n) = p0;                      % outside good price
    q0_c(n) = q0;                      % outside good quantity
    d_c(n)  = d;

    % Firm-specific values
    fprod_c(n,id(1:j)) = nextdi(2:j+1);             % productivity
    p_c(n,id(1:j)) = p(1:j);                        % price
    x_c(n,id(1:j)) = x(1:j);                        % rd investment
    q_c(n,id(1:j)) = q(1:j);                        % quantitiy sold
    rev_c(n,id(1:j)) = rev(1:j);
    pij_c(n,id(1:j)) = pr(1:j);                     % profits
    mgc_c(n,id(1:j)) = mgc(nextdi(2:j+1));          % marginal cost
    
    % New jump
    n = n + 1;
    if sh>0
        t = t - log(drw1(n))/sh;
    else
        % No more possible jumps (h=0), end loop.
        % This part will not run in general (either delta>0 or Q does not
        % have an absorbing state).
        t = T;
    end
end

% Compute discount rates.
temp1 = repmat(times_c(1:n-1),1,nf);
temp2 = tentry(1:n-1,1:nf);
discount = exp(-par.rho*(temp1 - temp2));

% Compile output.
out = struct();
out.times  = times_c(1:n-1,:);
out.discountf = discount;
out.fprod  = fprod_c(1:n-1,1:nf);
out.p      = p_c(1:n-1,1:nf);
out.x      = x_c(1:n-1,1:nf);
out.exit   = exit_c(1:n-1,:);
out.entry  = entry_c(1:n-1,:);
out.p0     = p0_c(1:n-1,:);
out.q0     = q0_c(1:n-1,:);
out.q      = q_c(1:n-1,1:nf);
out.rev    = rev_c(1:n-1,1:nf);
out.mgc    = mgc_c(1:n-1,1:nf);
out.pij    = pij_c(1:n-1,1:nf);
out.ecost  = ecost_c(1:n-1,:);
out.scrap  = scrap_c(1:n-1,:);
out.ecostj = ecostj_c(1:n-1,1:nf);
out.scrapj = scrapj_c(1:n-1,1:nf);
out.cs     = cs_c(1:n-1,:);
out.ps     = ps_c(1:n-1,:);
out.nf     = nf_c(1:n-1,:);
out.d      = d_c(1:n-1,:);
out.nfvec  = nf_c(1:n-1,:);
out.nf     = nf;
out.rdshock  = rdshock(1:n-1,1:nf);
out.lbdshock = lbdshock(1:n-1,1:nf);
out.n      = n-1;
out.flag   = flag;


end



%% Compute profits and quantities.
function [f,q,q0] = profits(par,w,p,p0,d)

% Note that here w is the index of productivity.
% Everything is computed for active firms only.

den = exp(-par.alpha*(p - p0));
s   = exp(-par.alpha*(p - p0))./(1 + sum(den)); % shares
q   = d*s;                                      % quantities
q0  = d*(1-sum(s));                             % outside good quantities
c   = par.mgc(w);                               % marginal cost

f = q.*(p - c);                                 % profits

end

%% Compute consumer surplus.
function cs = surplus(par,p,p0,d)

cs = d*(1/par.alpha)*log(sum(exp(-par.alpha*p))+exp(-par.alpha*p0));

end

%% Compute hazards of jumps.
function [hj,sh] = hazz(par,di,q,x,y)

% Computes 'hj', a vector of event-specific hazards, and 'sh', the hazard
% of a jump in the state of the industry due to any shock. 
% Note: hazard vectors are all column vectors.

in1 = (di(2:end)==par.M);
in2 = (di(2:end)==1);
in3 = (di(2:end)<par.M+1);
en  = find(di(2:end)==par.M+1,1);

h1 = (par.eta1*q)./(1+par.eta1*q);          % learning hazard
h1(in1) = 0;
h2 = par.alpha1*x(in3).^par.alpha2;         % r&d hazard
h2(in1)= 0;
h3 = ones(sum(in3),1)*par.delta;            % negative prod. shock hazard
h3(in2) = 0;
h4 = par.etax1*y(in3);                      % exit hazard
h5 = par.etae1*y(en);                       % entry hazard
h6 = par.gamma;                             % aggregate shock hazard 


if isempty(h1) 
    h1=zeros(in3,1);
end
if isempty(h2)
    h2=zeros(in3,1);
end
if isempty(h3)
    h3=zeros(in3,1);
end
if isempty(h4)
    h4=zeros(in3,1);
end
if isempty(h5)
    h5=0;
end

% h1 and h2 are not combined so that we can keep track of LBD vs R&D.
% Note the dim of hj is (4j+2,1).
hj = [h1;h2;h3;h4;h5;h6];
sh = sum(h1)+sum(h2)+sum(h3)+sum(h4)+h5+h6;

end

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
