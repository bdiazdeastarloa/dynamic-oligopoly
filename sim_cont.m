function [nf,times_c,ind_c,exit_c,entry_c,fprod_c,pij_c,price_c,p0_c,q0_c,mgc_c,quant_c,reve_c,ecostj_c,scrapj_c,ecost_c,scrap_c,csurplus_c,aggprof_c,nfirms_c,flag] = sim_cont(par,T,P,E,X,wstart,w0start,sim)

% =========================================================================
% Function: sim_cont
%
% Simulates arrival times and state jumps of the continuous time model
% during [0,T];
% Inputs are: the parameter structure, the number of periods (T), policy 
% functions (prices, entry, exit), and the initial industry structure (w0). 
% 
% Written by BDA @ PSU July 2014
% =========================================================================

draws1 = par.randset1(:,sim);
draws2 = par.randset2(:,sim);
draws3 = par.randset3(:,sim);
I      = eye(par.p0_size);            % identity for transition prob
Q      = par.Q/par.gamma;             % p0 transition matrix

times_c = zeros(T*20,1);             % arrival times
ind_c   = zeros(T*20,par.nfirms+1);  % industry structures        
exit_c  = zeros(T*20,1);             % exit events
entry_c = zeros(T*20,1);             % entry events
fprod_c  = zeros(T*20,20);           % firms' productivities
price_c = zeros(T*20,20);            % firms' prices
p0_c    = zeros(T*20,1);             % outside good price
q0_c    = zeros(T*20,1);             % outside good quantity
quant_c = zeros(T*20,20);            % firms' quantities
reve_c  = zeros(T*20,20);            % firms' revenues
mgc_c   = zeros(T*20,20);            % firms' marginal cost
pij_c   = zeros(T*20,20);            % firms' profits
ecostj_c = zeros(T*20,20);        
scrapj_c = zeros(T*20,20);
ecost_c  = zeros(T*20,1);            % entry costs
scrap_c  = zeros(T*20,1);            % exit scrap values
csurplus_c = zeros(T*20,1);          % consumer surplus
aggprof_c  = zeros(T*20,1);          % aggregate profits
nfirms_c = zeros(T*20,1);            % number of active firms

locw  = wstart;                      % initial industry structure            
w0 = w0start;
in = (locw>0);                       % index of active firms
j  = sum(in);
nf = j;
w  = qencode(par,locw);              % encoded ind structure
w  = w + par.dostates*(w0-1);        % state position
p0 = par.p0vec(w0);

price = P(w,in);
ex    = X(w,in);
en    = E(w,:);
[h1,h2,h3,h4,h] = hazz(par,locw,price,p0,ex,en);
[Pi,q,q0] = profits(par,locw(in),price',p0);
csurp   = surplus(par,price,p0);
rev     = price'.*q;
aggprof = sum(Pi);
n = 1;                                % first "jump" is initial state
nextw = locw;

% Record initial state
ind_c(n,:)    = [w0 nextw'];          % industry structure
csurplus_c(n) = csurp;                % consumer surplus
aggprof_c(n)  = aggprof;              % aggregate profits
nfirms_c(n)   = j;                    % # of active firms
p0_c(n) = p0;
q0_c(n) = q0;
id = zeros(par.nfirms,1);             % firms' positions
id(1:j) = (1:1:j);

% Firm-specific values
for ff=1:j
    fprod_c(n,id(ff))  = nextw(ff);                % productivity
    price_c(n,id(ff)) = price(ff);                 % price
    quant_c(n,id(ff)) = q(ff);                     % quantitiy sold
    reve_c(n,id(ff))  = rev(ff);
    pij_c(n,id(ff))   = Pi(ff);                    % profits
    mgc_c(n,id(ff))   = par.mgc(nextw(ff));        % marginal cost
end

flag = 0; 
n = n + 1;
t = - log(draws1(n))/h;                            % first actual jump
if t==Inf
    disp('WARNING: t=Inf - no first jump. Something is wrong.');
    flag = 1;
end
   
while t<T
    % There's a jump!
    times_c(n) = t;                   
    % Where to jump?
    % (replicates 'mnrnd(1,x)' to get the random draw for entry/exit)
    x    = [h1,h2,h3,h4]./h;                       % jump probabilities
    prob = draws2(n);
    edge = [zeros(1,1) cumsum(x,2)];
    draw = histc(prob,edge);
    draw(:,end) = [];
    draw = find(draw==1);
    % Outside good changed
    if draw==size(x,2)
        % (replicates 'mnrnd(1,x)')
        trans = I(w0,:)*Q;                 % p0 conditional transition prob
        prob = draws3(n);
    	temp = [zeros(1,1) cumsum(trans,2)];
        draw = histc(prob,temp);
        draw(:,end) = [];
        w0 = find(draw==1);
        p0 = par.p0vec(w0);
    % Entry
    elseif draw==size(x,2)-1                     
        ii = j+1;                          % new firm index
        wentry = max(nextw(1) - 2,1);      % entry draw
        nextw(ii) = wentry;
        nf = nf + 1;                       % add one firm to total firms
        id(ii) = nf;                       % new firm ID
        % Update entry event and entry cost
        entry_c(n) = 1;
        ecost_c(n) = par.entrylo + prob*(par.entryhi-par.entrylo);
        ecostj_c(n,id(ii)) = ecost_c(n);
    % Exit
    elseif draw>j                          % 1st j cells are omega jumps         
        ii = draw-j;                       % which firm exits
        nextw(ii) = 0; 
        % Update exit event and exit cost
        exit_c(n)  = 1;
        scrap_c(n) = par.scraplo + prob*(par.scraphi-par.scraplo);
        scrapj_c(n,id(ii)) = scrap_c(n);
        id(ii) = 0;                        % free ID slot
    % Productivity change
    else                                 
        ii = draw;                         % which firm goes up
        nextw(ii) = nextw(ii)+1;
    end
    
    temp  = sortrows(flipud([nextw,id]),-1);  % sort decreasing  
    nextw = temp(:,1);
    id    = temp(:,2);              
       
    % Update everything
    in = (nextw>0); 
    j  = sum(in);
    if j==0
        disp('WARNING: empty industry.');
        flag=1;
        break
    end
    w  = qencode(par,nextw);
    w  = w + par.dostates*(w0-1);
    price = P(w,in);
    ex  = X(w,in);
    en  = E(w,:);
    [h1,h2,h3,h4,h] = hazz(par,nextw,price,p0,ex(in),en);
    [Pi,q,q0] = profits(par,nextw(in),price',p0);
    csurp  = surplus(par,price,p0);
    rev    = price'.*q;
    aggprof = sum(Pi);
    
    % Record outcomes
    ind_c(n,:) = [w0 nextw'];                   % industry structure
    csurplus_c(n) = csurp;                      % consumer surplus
    aggprof_c(n)  = aggprof;                    % aggregate profits
    nfirms_c(n) = j;                            % # of active firms
    p0_c(n) = p0;                               % outside good price
    q0_c(n) = q0;                               % outside good quantity
    
    % Firm-specific values
    for ff=1:j
        fprod_c(n,id(ff)) = nextw(ff);                 % productivity
        price_c(n,id(ff)) = price(ff);                 % price
        quant_c(n,id(ff)) = q(ff);                     % quantitiy sold
        reve_c(n,id(ff))  = rev(ff);                   % revenues
        pij_c(n,id(ff))   = Pi(ff);                    % profits
        mgc_c(n,id(ff))   = par.mgc(nextw(ff));        % marginal cost
    end
    
    % New jump
    n = n + 1;
    if h>0
        t = t - log(draws1(n))/h;
    else
        % No more possible jumps (h=0), end loop 
        % This part should not run if p0 does not have absorbing state
        % because there will always be a jump
        t = T;
    end
end

times_c  = times_c(1:n-1,:);
ind_c    = ind_c(1:n-1,:);
exit_c   = exit_c(1:n-1,:);
entry_c  = entry_c(1:n-1,:);
fprod_c  = fprod_c(1:n-1,1:nf);
price_c  = price_c(1:n-1,1:nf);
p0_c     = p0_c(1:n-1,:);
q0_c     = q0_c(1:n-1,:);
quant_c  = quant_c(1:n-1,1:nf);
reve_c   = reve_c(1:n-1,1:nf);
mgc_c    = mgc_c(1:n-1,1:nf);
pij_c    = pij_c(1:n-1,1:nf);
ecost_c  = ecost_c(1:n-1,:);
scrap_c  = scrap_c(1:n-1,:);
ecostj_c = ecostj_c(1:n-1,1:nf);
scrapj_c = scrapj_c(1:n-1,1:nf);
csurplus_c = csurplus_c(1:n-1,:);
aggprof_c  = aggprof_c(1:n-1,:);
nfirms_c   = nfirms_c(1:n-1,:);
