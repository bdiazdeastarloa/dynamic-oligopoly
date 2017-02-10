function [aggdata,paneldata] = sim_disc(par,outc,T)

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
% 
% sim_disc: fits simulation output of 'sim_cont' into a discrete time grid.
%
% Output is:
% - aggdata cell array: industry (aggregate) stats.
% - paneldata cell array: firm-specific stats (unbalanced panel).
%
% Written by Bernardo Diaz de Astarloa @ PSU 2014
% =========================================================================

if outc.flag==1
    aggdata = ones(T,14);
    paneldata = ones(T*outc.nf,20);
    ff        = repmat((1:outc.nf)',T,1);
    paneldata = [ff paneldata];
else
    % Load output from simulation.
    times_c = outc.times;
    discountf = outc.discountf;
    fprod_c = outc.fprod;
    x_c     = outc.x;
    exit_c  = outc.exit;
    entry_c = outc.entry;
    p0_c    = outc.p0;
    q0_c    = outc.q0;
    q_c     = outc.q;
    rev_c   = outc.rev;
    mgc_c   = outc.mgc;
    pij_c   = outc.pij.*discountf;
    ecost_c = outc.ecost;
    scrap_c = outc.scrap;
    ecostj_c= outc.ecostj;
    scrapj_c= outc.scrapj;
    cs_c    = outc.cs;
    ps_c    = outc.ps;
    psx_c   = outc.psx;
    d_c     = outc.d;
    nfvec   = outc.nfvec;
    rd_c    = outc.rdshock;
    lbd_c   = outc.lbdshock;
    nf      = outc.nf;
    
    ff     = zeros(T,1);        % # of firms
    exit   = zeros(T,1);        % # of exiting firms
    entry  = zeros(T,1);        % # of entrants
    p0     = zeros(T,1);        % imports price
    q0     = zeros(T,1);        % imported quantities
    sh0    = zeros(T,1);        % imports share
    d      = zeros(T,1);        % demand
    cs     = zeros(T,1);        % consumer surplus
    ps     = zeros(T,1);        % aggregate profits
    psx    = zeros(T,1);        % aggregate profits net of R&D
    scrap  = zeros(T,1);        % total scrap value
    ecost  = zeros(T,1);        % total entry costs

    w      = zeros(T,nf);        % productivity 
    q      = zeros(T,nf);        % quantities
    rev    = zeros(T,nf);        % revenues
    p      = zeros(T,nf);        % average price (rev/q)
    x      = zeros(T,nf);        % r&d investment
    pij    = zeros(T,nf);        % discounted profits
    mgc    = zeros(T,nf);        % marginal cost
    scrapj = zeros(T,nf);        % exit (firm-specific)
    ecostj = zeros(T,nf);        % entry (firm-specific)
    sh     = zeros(T,nf);        % market share
    rd     = zeros(T,nf);        % r&d steps
    lbd    = zeros(T,nf);        % lbd steps

    paneldata = [];
    FT        = cell(1,nf);
    for j=1:nf 
        FT{j} = zeros(T,21);     % Create holder for each firm
        FT{j}(:,1) = j;          % Firm id       
        FT{j}(:,2) = (1:1:T);    % Time period
    end

    % Replace productivity indices with their values.
    fprod_c((fprod_c==0))= par.M+1;
    temp = [par.w 0]';
    fprod_c = temp(fprod_c);


    %% Get aggregate and firm trajectories.

    t_lag = 0;
    for t = 1:T
        t_ind = find(times_c(:,1)<t,1,'last');
        if isempty(t_ind)
            t_ind = t_lag;
        end
        if t_ind == t_lag
            if t==1
                % No jumps in t=1: keep initial state.
                ff(t) = nfvec(1);
                p0(t) = p0_c(1);
                q0(t) = q0_c(1);
                d(t)  = d_c(1);
                sh0(t) = q0_c(1)/d_c(1);
                cs(t) = cs_c(1);
                ps(t) = ps_c(1);
                psx(t) = psx_c(1);

                w(t,:)   = par.w(fprod_c(1,:));          
                q(t,:)   = q_c(1,:);
                rev(t,:) = rev_c(1,:);
                p(t,:)   = rev(t)./q(t);
                x(t,:)   = x_c(1,:);
                pij(t,:) = pij_c(1,:); 
                mgc(t,:) = mgc_c(1,:);  
                ecostj(t,:) = ecostj_c(1,:);
                scrapj(t,:) = scrapj_c(1,:);
                sh(t,:)  = q(t)./d(t);
            else
                % No jumps in t>1: same as last period.
                ff(t) = ff(t-1);
                p0(t) = p0(t-1);
                q0(t) = q0(t-1);
                d(t)  = d(t-1);
                sh0(t) = sh0(t-1);
                cs(t) = cs(t-1);
                ps(t) = ps(t-1);
                psx(t) = psx(t-1);

                % Adjust if there was exit in last jump, last period.
                ex_in = (scrapj(t-1,1:nf)>0);   % which firms exited: they have zero everything.           
                w(t,ex_in)   = 0;
                q(t,ex_in)   = 0;        
                rev(t,ex_in) = 0;
                p(t,ex_in)   = 0;
                x(t,ex_in)   = 0;
                pij(t,ex_in) = 0;
                mgc(t,ex_in) = 0;
                sh(t,ex_in)  = 0;

                cont_in = (ex_in==0);           % which firms continued: same as last.
                w(t,cont_in)   = w(t-1,cont_in);
                q(t,cont_in)   = q(t-1,cont_in);        
                rev(t,cont_in) = rev(t-1,cont_in);
                p(t,cont_in)   = p(t-1,cont_in);
                x(t,cont_in)   = x(t-1,cont_in);
                pij(t,cont_in) = pij(t-1,cont_in).*exp(-par.rho);
                mgc(t,cont_in) = mgc(t-1,cont_in);
                sh(t,cont_in)  = sh(t-1,cont_in);

                % Adjust if there was entry in last jump, last period.
                % (Do not duplicate entry cost.)
                e_in = (ecostj(t-1,1:nf)>0);                   
                ecostj(t,e_in) = 0;                   

            end % if t==1
        else 
            % There was a jump.        
            % Weighted (integrated) variables.
            if t==1      % t_lag=0 for t=1, cannot index with 0.
                weight = [times_c(t_lag+2:t_ind);t] - times_c(t_lag+1:t_ind);
                p0(t)  = weight'*p0_c(t_lag+1:t_ind);
                q0(t)  = weight'*q0_c(t_lag+1:t_ind);
                cs(t)  = weight'*cs_c(t_lag+1:t_ind);
                ps(t)  = weight'*ps_c(t_lag+1:t_ind);
                psx(t) = weight'*psx_c(t_lag+1:t_ind);
                d(t)   = weight'*d_c(t_lag+1:t_ind);
                sh0(t) = q0(t)/d(t);

                q(t,:)   = weight'*q_c(t_lag+1:t_ind,:); 
                rev(t,:) = weight'*rev_c(t_lag+1:t_ind,:);
                x(t,:)   = weight'*x_c(t_lag+1:t_ind,:);
                pij(t,:) = weight'*pij_c(t_lag+1:t_ind,:);

                for j = 1:nf
                    active = (fprod_c(t_lag+1:t_ind,j)>0);
                    if sum(active)==0 
                        w(t,j)=0;
                        mgc(t,j)=0;
                    else
                        temp = fprod_c(t_lag+1:t_ind,j);
                        w(t,j) = (weight(active)'*temp(active))/sum(weight(active));
                        temp = mgc_c(t_lag+1:t_ind,j);
                        mgc(t,j) = (weight(active)'*temp(active))/sum(weight(active));
                    end
                end
            else
                weight = [times_c(t_lag+1:t_ind);t] - [t-1;times_c(t_lag+1:t_ind)];
                p0(t)  = weight'*p0_c(t_lag:t_ind);
                q0(t)  = weight'*q0_c(t_lag:t_ind);
                cs(t)  = weight'*cs_c(t_lag:t_ind);
                ps(t)  = weight'*ps_c(t_lag:t_ind);
                psx(t) = weight'*psx_c(t_lag:t_ind);
                d(t)   = weight'*d_c(t_lag:t_ind);
                sh0(t) = q0(t)/d(t);

                q(t,:)   = weight'*q_c(t_lag:t_ind,:); 
                rev(t,:) = weight'*rev_c(t_lag:t_ind,:);
                x(t,:)   = weight'*x_c(t_lag:t_ind,:);
                pij(t,:) = weight'*pij_c(t_lag:t_ind,:);
                
                for j = 1:nf
                    active = (fprod_c(t_lag:t_ind,j)>0);
                    if sum(active)==0
                        w(t,j)=0;
                        mgc(t,j)=0;
                    else
                        temp = fprod_c(t_lag:t_ind,j);
                        w(t,j) = (weight(active)'*temp(active))/sum(weight(active));
                        temp = mgc_c(t_lag:t_ind,j);
                        mgc(t,j) = (weight(active)'*temp(active))/sum(weight(active));
                    end
                end
            end
            
            % Non-weighted variables.
            % Start from t_lag since the period starts with last period's
            % last jump # of firms.
            if t==1
                ff(t)    = max(nfvec(t_lag+1:t_ind));
            else
                ff(t)    = max(nfvec(t_lag:t_ind));
            end
            exit(t)  = sum(exit_c(t_lag+1:t_ind));
            entry(t) = sum(entry_c(t_lag+1:t_ind));
            scrap(t) = sum(scrap_c(t_lag+1:t_ind));
            ecost(t) = sum(ecost_c(t_lag+1:t_ind));

            sh(t,:) = q(t,:)./d(t);
            p(t,:)  = rev(t,:)./q(t,:);
            ecostj(t,:) = sum(ecostj_c(t_lag+1:t_ind,:));
            scrapj(t,:) = sum(scrapj_c(t_lag+1:t_ind,:));
            rd(t,:) = sum(rd_c(t_lag+1:t_ind,:));
            lbd(t,:) = sum(lbd_c(t_lag+1:t_ind,:));
            
            % Update starting point for time period tracker.
            t_lag = t_ind;

            % Store firm-level data in cells.
            for j = 1:nf         
                FT{j}(t,3)  = w(t,j);                   % productivity
                FT{j}(t,4)  = p(t,j);                   % firm price
                FT{j}(t,5)  = x(t,j);                   % R&D investment
                FT{j}(t,6)  = q(t,j);                   % firm quantities
                FT{j}(t,7)  = rev(t,j);                 % firm revenues
                FT{j}(t,8)  = pij(t,j);                 % firm profits
                FT{j}(t,9)  = mgc(t,j);                 % marginal cost
                FT{j}(t,10) = sh(t,j);                  % firm market share
                FT{j}(t,11) = ecostj(t,j) ;             % entry cost
                FT{j}(t,12) = scrapj(t,j);              % scrap value
                FT{j}(t,13) = rd(t,j);                  % R&D jumps
                FT{j}(t,14) = lbd(t,j);                 % LBD jumps
            end
        end % if t_ind == t_lag
    end % for t=1:T


    %% Construct panel data.

    for j=1:nf
        active = (FT{j}(:,3)>0);                        % active periods (positive productivity) 
        FT{j}  = FT{j}(active,:);                       % keep active years only

        FT{j}(:,15) = [NaN(1); FT{j}(1:end-1,4)];       % lagged prices
        FT{j}(:,16) = [NaN(1); FT{j}(1:end-1,6)];       % lagged shipments
        FT{j}(:,17) = [NaN(1); FT{j}(1:end-1,7)];       % lagged revenues
        FT{j}(:,18) = [NaN(1); FT{j}(1:end-1,10)];      % lagged shares 

        temp = 1:1:T;
        FT{j}(:,19) = p0(temp(active));                 % foreign price
        FT{j}(:,20) = [NaN(1); FT{j}(1:end-1,19)];      % lagged foreign price
        FT{j}(:,21) = [NaN(1); FT{j}(1:end-1,5)];       % lagged R&D
        
        FT{j}(:,22) = cumsum(FT{j}(:,6));               % cumulative shipments
        FT{j}(:,23) = cumsum(FT{j}(:,5));               % cumulative r&d
        
        if size(FT{j}(:,1),1)<2
            FT{j}(:,22) = NaN;
        else
            FT{j}(:,22) = [NaN(2,1); FT{j}(1:end-2,6)];     % lagged 2 shipments
        end

        % Panel data with simulated variables.
        paneldata = [paneldata; FT{j}];                  %#ok
    end

    %% Construct aggregate data.

    % Selected firm outcomes.
    w       = zeros(T,1);
    avprice = zeros(T,1);
    c2      = zeros(T,1);
    avrd    = zeros(T,1);
    totrd   = zeros(T,1);
    avmgc   = zeros(T,1);
    rdjump  = zeros(T,1);
    lbdjump = zeros(T,1);
    
    if nf>0
        w_temp       = sortrows([paneldata(:,2),paneldata(:,3)],1);
        avprice_temp = sortrows([paneldata(:,2),paneldata(:,4)],1);          
        c2_temp      = sortrows([paneldata(:,2),paneldata(:,10)],[1 -2]);     % sort by decreasing share
        x_temp       = sortrows([paneldata(:,2),paneldata(:,5)],1);
        mgc_temp     = sortrows([paneldata(:,2),paneldata(:,9)],1);
        rdjump_temp  = sortrows([paneldata(:,2),paneldata(:,13)],1);
        lbdjump_temp = sortrows([paneldata(:,2),paneldata(:,14)],1);

        for t = 1:T
            % Average productivity (unweighted).
            w(t) = mean(w_temp(w_temp(:,1)==t,2));     
            % Average prices.
            avprice(t) = mean(avprice_temp(avprice_temp(:,1)==t,2));   
            % Concentration (C2 index).
            temp = c2_temp(c2_temp(:,1)==t,2);            
            c2(t) = sum(temp(1:min(size(temp,1),2))); 
            % Average and total R&D.
            avrd(t) = mean(x_temp(x_temp(:,1)==t,2));
            totrd(t) = sum(x_temp(x_temp(:,1)==t));
            % Average marginal cost.
            avmgc(t) = mean(mgc_temp(mgc_temp(:,1)==t,2));
            % # of RD and LBD jumps.
            rdjump(t) = mean(rdjump_temp(rdjump_temp(:,1)==t,2));
            lbdjump(t) = mean(lbdjump_temp(lbdjump_temp(:,1)==t,2));
        end
    end

    % Aggregate simulated variables.
    aggdata = [ff,entry,exit,cs,ps,psx,scrap,ecost,p0,q0,sh0,d,w,avprice,c2,avrd,totrd,avmgc,rdjump,lbdjump];
end
