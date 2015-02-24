function [aggdata,paneldata,averages] = sim_disc(par,nf,times_c,exit_c,entry_c,fprod_c,pij_c,p0_c,q0_c,mgc_c,quant_c,rev_c,ecostj_c,scrapj_c,ecost_c,scrap_c,csurplus_c,aggprof_c,nfirms_c,T)

% =========================================================================
% Function: sim_disc
%
% Fits simulation output of 'sim_cont' into a discrete time grid.
%
% Output is:
% - aggdata cell array: industry (aggregate) stats.
% - paneldata cell array: firm-specific stats (unbalanced panel).
%
% Written by BDA @ PSU July 2014
% =========================================================================

nfirms   = zeros(T,1);      % # of firms
exit     = zeros(T,1);      % # of exiting firms
entry    = zeros(T,1);      % # of entrants
p0       = zeros(T,1);      % imports price
q0       = zeros(T,1);      % imported quantities
share0   = zeros(T,1);      % imports share
csurp    = zeros(T,1);      % consumer surplus
aggpi    = zeros(T,1);      % aggregate profits
scrap    = zeros(T,1);      % total scrap value
ecost    = zeros(T,1);      % total entry costs

omega  = zeros(T,1);        % productivity 
q      = zeros(T,1);        % quantities
rev    = zeros(T,1);        % revenues
p      = zeros(T,1);        % average price (rev/q)
pij    = zeros(T,1);        % profits
mgc    = zeros(T,1);        % marginal cost
scrapj = zeros(T,1);        % exit (firm-specific)
ecostj = zeros(T,1);        % entry (firm-specific)
share  = zeros(T,1);        % market share

M = par.M;
w = par.wvec';
paneldata = [];
FT        = cell(1,nf);

        
%% Get aggregate industry statistics

t_lag = 0;
for t = 1:T
    t_ind = find(times_c(:,1)<t,1,'last');
    if isempty(t_ind)
        t_ind=t_lag;
    end
    if t_ind == t_lag
        if t==1
            % No jumps in period 1: keep initial state
            nfirms(t) = nfirms_c(1);
            p0(t)     = p0_c(1);
            q0(t)     = q0_c(1);
            share0(t) = q0(t)/M;
            csurp(t)  = csurplus_c(1);
            aggpi(t)  = aggprof_c(1);
        else
            % Periods>1 with no jumps: same as last period
            nfirms(t) = nfirms(t-1);
            p0(t)     = p0(t-1);
            q0(t)     = q0(t-1);
            share0(t) = q0(t)/M;
            csurp(t)  = csurp(t-1);
            aggpi(t)  = aggpi(t-1);
        end 
    else 
        % There was a jump        
        % Weighted (integrated) variables
        if t==1      % t_lag=0 for t=1, cannot index with 0
            weight   = [times_c(t_lag+2:t_ind);t] - times_c(t_lag+1:t_ind);
            p0(t)    = weight'*p0_c(t_lag+1:t_ind);
            q0(t)    = weight'*q0_c(t_lag+1:t_ind);
            csurp(t) = weight'*csurplus_c(t_lag+1:t_ind);
            aggpi(t) = weight'*aggprof_c(t_lag+1:t_ind);
        else
            weight   = [times_c(t_lag+1:t_ind);t] - [t-1;times_c(t_lag+1:t_ind)];
            p0(t)    = weight'*[p0_c(t_lag);p0_c(t_lag+1:t_ind)];
            q0(t)    = weight'*[q0_c(t_lag);q0_c(t_lag+1:t_ind)];
            csurp(t) = weight'*[csurplus_c(t_lag);csurplus_c(t_lag+1:t_ind)];
            aggpi(t) = weight'*[aggprof_c(t_lag);aggprof_c(t_lag+1:t_ind)];
        end
        % Non-weighted variables   
        nfirms(t) = max(nfirms_c(t_lag+1:t_ind));
        exit(t)   = sum(exit_c(t_lag+1:t_ind));
        entry(t)  = sum(entry_c(t_lag+1:t_ind));
        scrap(t)  = sum(scrap_c(t_lag+1:t_ind));
        ecost(t)  = sum(ecost_c(t_lag+1:t_ind));
        share0(t) = q0(t)/M;
        t_lag = t_ind;
    end
end
 
%% Get firm trajectories
%  This is to be changed. Loop should be for t, then for j. 
%  Jumps are industry-wise, so doesn't make sense to have the j loop first,
%  it makes it too redundant.
for j = 1:nf
    FT{j} = zeros(T,17);            % Create on panel holder for each firm
    FT{j}(:,1) = j;                 % Firm ID
    FT{j}(:,2) = (1:1:T);           % Period index
    t_lag = 0;
    for t = 1:T        
        t_ind = find(times_c(:,1)<t,1,'last');
        if isempty(t_ind)
            t_ind=t_lag;
        end
        % No jumps
        if t_ind == t_lag
            if t==1
                % No jumps in period 1: keep initial state
                omega(t)  = w(fprod_c(1,j));          
                q(t)      = quant_c(1,j);
                rev(t)    = rev_c(1,j);
                p(t)      = rev(t)/q(t); 
                pij(t)    = pij_c(1,j); 
                mgc(t)    = mgc_c(1,j);  
                ecostj    = ecostj_c(1,j);
                scrapj(t) = scrapj_c(1,j);
                share(t)  = q(t)/M;
            else
                % Periods>1 with no jumps: same as last period
                if scrapj(t-1)>0                        % exit in last jump
                    omega(t) = 0;
                    q(t)     = 0;        
                    rev(t)   = 0;
                    p(t)     = 0;        
                    pij(t)   = 0;                 
                    mgc(t)   = 0;
                    share(t) = 0;
                else                                
                    omega(t) = omega(t-1);
                    q(t)     = q(t-1);        
                    rev(t)   = rev(t-1);
                    p(t)     = p(t-1);
                    pij(t)   = pij(t-1);                 
                    mgc(t)   = mgc(t-1);
                    share(t) = share(t-1);
                end
                if ecostj(t-1)>0                   % entry in last jump
                    ecostj(t) = 0;                 % do not repeat entry
                end
            end    
        else
            % There was a jump            
            % Weighted variables
            if t==1      % t_lag=0 for t=1, cannot index with 0
                weight   = [times_c(t_lag+2:t_ind);t] - times_c(t_lag+1:t_ind);
                q(t)     = weight'*quant_c(t_lag+1:t_ind,j); 
                rev(t)   = weight'*rev_c(t_lag+1:t_ind,j);
                pij(t)   = weight'*pij_c(t_lag+1:t_ind,j);
                mgc(t)   = weight'*mgc_c(t_lag+1:t_ind,j);
                w_temp   = fprod_c(t_lag+1:t_ind,j);
                active   = find(w_temp>0);
            else
                weight   = [times_c(t_lag+1:t_ind);t] - [t-1;times_c(t_lag+1:t_ind)];
                q(t)     = weight'*quant_c(t_lag:t_ind,j); 
                rev(t)   = weight'*rev_c(t_lag:t_ind,j);
                pij(t)   = weight'*pij_c(t_lag:t_ind,j);
                mgc(t)   = weight'*mgc_c(t_lag:t_ind,j);
                w_temp   = fprod_c(t_lag:t_ind,j);
                active   = find(w_temp>0);
            end
            if isempty(active)==1
                omega(t) = 0;
            else
                w_temp(active) = w(w_temp(active));
                omega(t) = (weight'*w_temp)/sum(weight(active));
            end
            
            % Unweighted variables
            share(t) = q(t)/M;
            p(t)     = rev(t)/q(t);

            % Entry cost
            ecostj(t) = sum(ecostj_c(t_lag+1:t_ind,j));
            
            % Exit scrap value
            scrapj(t) = sum(scrapj_c(t_lag+1:t_ind,j));
        end           
        t_lag = t_ind;
        
        FT{j}(t,3)  = omega(t);
        FT{j}(t,4)  = p(t);
        FT{j}(t,5)  = q(t);
        FT{j}(t,6)  = rev(t);
        FT{j}(t,7)  = pij(t);
        FT{j}(t,8)  = mgc(t);
        FT{j}(t,9)  = share(t);
        FT{j}(t,10) = ecostj(t) ;
        FT{j}(t,11) = scrapj(t);
    end

    active = (FT{j}(:,3)>0);                        % active periods (positive productivity) 
    FT{j}  = FT{j}(active,:);                       % keep active years only
    
    FT{j}(:,12) = cumsum(FT{j}(:,6));               % cumulative shipments
    FT{j}(:,13) = [NaN(1); FT{j}(1:end-1,4)];       % lagged prices
    FT{j}(:,14) = [NaN(1); FT{j}(1:end-1,5)];       % lagged shipments
    FT{j}(:,15) = [NaN(1); FT{j}(1:end-1,6)];       % lagged revenues
    FT{j}(:,16) = [NaN(1); FT{j}(1:end-1,8)];       % lagged shares 
   
    if sum(active)>2
        FT{j}(:,17) = [NaN(2,1); FT{j}(1:end-2,6)]; % lagged*2 revenues
        FT{j}(:,18) = [NaN(2,1); FT{j}(1:end-2,8)]; % lagged*2 shares
        
        temp = zeros(size(FT{j}(:,1),1),1);
        for i=1:size(FT{j}(:,1),1)-2;
            temp(i+2) = sum(FT{j}(i:i+2,6));
        end
        FT{j}(:,18) = [NaN(2,1); temp(3:end,1)];    % "t:t-2" cum shipments
    else
        FT{j}(:,16) = NaN(sum(active),1);           % lagged*2 revenues
        FT{j}(:,17) = NaN(sum(active),1);           % lagged*2 shares
        FT{j}(:,18) = NaN(sum(active),1);           % "t:t-2" cum shipments
    end
    
    % Compile firm-level simulated variables
    paneldata = [paneldata; FT{j}];                  %#ok
end

% Compile aggregate simulated variables
aggdata = [nfirms,entry,exit,csurp,aggpi,scrap,ecost,p0,q0,share0];

% Average of some firm outcomes (conevenient to do it here and allow
% paralellization of simulations)
avprice  = zeros(T,1);
c2       = zeros(T,1);
omega    = zeros(T,1);

omega_temp   = sortrows([paneldata(:,2),paneldata(:,3)],1);
avprice_temp = sortrows([paneldata(:,2),paneldata(:,4)],1);          
c2_temp      = sortrows([paneldata(:,2),paneldata(:,8)],[1 -2]);     % sort by decreasing share
for t = 1:T
    % Average productivity (unweighted)
    omega(t) = mean(omega_temp(omega_temp(:,1)==t,2));     
    % Average prices
    avprice(t) = mean(avprice_temp(avprice_temp(:,1)==t,2));   
    % Concentration (C2 index)
    temp = c2_temp(c2_temp(:,1)==t,2);            
    c2(t) = sum(temp(1:min(size(temp,1),2)));    
end

averages = [omega,avprice,c2];

