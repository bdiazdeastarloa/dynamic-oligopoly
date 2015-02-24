function MM = moms(par,P,X,E,solveflag)

% =========================================================================
%
% Creates moments to match of the continuous time learning by doing model.
% Calls:
% - sim_cont: simulates arrival times and outcomes.
% - sim_disc: fits simulation in a discrete time grid. 
%
% Written by BDA @ PSU July 2014.
% =========================================================================

%% Simulation settings
T = par.T;              % # of periods to simulate
S = par.S;              % # of simulations to run

wstart  = par.wstart;   % initial industry structure
w0start = par.w0start;  % starting outside good price

learning  = zeros(1,S);
share_ar1 = zeros(1,S);
reve_ar1  = zeros(1,S);
exitrate  = zeros(1,S);
entryrate = zeros(1,S);
avp       = zeros(1,S);

if solveflag==1
    learning  = ones(1,S)*100;
    share_ar1 = ones(1,S)*100;
    reve_ar1  = ones(1,S)*100;
    exitrate  = ones(1,S)*100;
    entryrate = ones(1,S)*100;
    avp       = ones(1,S)*100;
else   
    for sim=1:S
        % Simulate industry: arrival times and jumps
        [nf,times_c,~,exit_c,entry_c,fprod_c,pij_c,~,p0_c,q0_c,mgn_c,quant_c,rev_c,ecostj_c,scrapj_c,ecost_c,scrap_c,csurplus_c,aggprof_c,nfirms_c,simflag] = sim_cont(par,T,P,E,X,wstart,w0start,sim);
        
        if simflag==1
            disp('Something went wrong in simulation');
            learning  = ones(1,S)*100;
            share_ar1 = ones(1,S)*100;
            reve_ar1  = ones(1,S)*100;
            exitrate  = ones(1,S)*100;
            entryrate = ones(1,S)*100;
            avp       = ones(1,S)*100;
            break
        else   
            % Fit to discrete time grid (years)
            % AD: aggregate data | PD: panel data
            [agg,panel,~] = sim_disc(par,nf,times_c,exit_c,entry_c,fprod_c,pij_c,p0_c,q0_c,mgn_c,quant_c,rev_c,ecostj_c,scrapj_c,ecost_c,scrap_c,csurplus_c,aggprof_c,nfirms_c,T);

            % Define variables for moments
            % Simulated data are sorted by id then period
            id     = panel(:,1);
            year   = panel(:,2);
            price  = panel(:,4);
            rev    = panel(:,6);
            share  = panel(:,8);

            lagrev   = panel(:,14);
            lagshare = panel(:,15);
            
            cumq1    = panel(:,11);
            cumq2    = panel(:,18);
            
            entry = agg(:,2);
            exit  = agg(:,3);
            
            const = ones(size(panel,1),1);

            % ---------------------------------------------------------------------
            % Compute simulated moments
            % ---------------------------------------------------------------------
            
            % ---------------------------------------------------------------------
            % Leaning curve
            ii = (cumq1>0);
            %ii = (cumq2>0);
            y  = log(price(ii,:));
            
            % Generate firm fixed effects matrix
            temp   = eye(nf);
            fix_id = zeros(size(y,1),nf); 
            tempid = id(ii);
            for i=1:size(y,1)
                fix_id(i,:) = (temp(tempid(i),:));
            end
            % Check for columns of zeros
            check  = (sum(fix_id)>0);
            fix_id = fix_id(:,check); 
            x = [const(ii) log(cumq1(ii)) fix_id(:,2:end)];
            % x = [const(ii) log(cumq2(ii)) fix_id(:,2:end)];
            
            [b,~,~] = regress(y,x);
            % rMSE = sqrt(sum(r.^2)/size(y,1));
            learning(:,sim) = b(2);
            % ---------------------------------------------------------------------
            
            % ---------------------------------------------------------------------
            % Exit rate
            exitrate(:,sim) = mean(exit);
            % Entry rate
            entryrate(:,sim) = mean(entry);
            % ---------------------------------------------------------------------
            
            % ---------------------------------------------------------------------
            % Log shares AR1
            ii = (lagshare>0);
            y  = log(share(ii,:));
            x  = [const(ii) year(ii) log(lagshare(ii))];
            [b,~,~] = regress(y,x);
            % rMSE    = sqrt(sum(r.^2)/(size(y,1)-2));
            share_ar1(:,sim) = b(3);
            % ---------------------------------------------------------------------
            
            % ---------------------------------------------------------------------
            % Log revenues AR1
            ii = (lagrev>0);
            y  = log(rev(ii,:));
            x  = [const(ii) year(ii) log(lagrev(ii))];
            [b,~,~] = regress(y,x);
            % rMSE    = sqrt(sum(r.^2)/(size(y,1)-2));
            reve_ar1(:,sim) = b(3);
            % ---------------------------------------------------------------------
            
            % ---------------------------------------------------------------------
            % Mean price
            avp(:,sim) = mean(price(price>0));
            % ---------------------------------------------------------------------
        end
    end
end

learning  = mean(learning,2);
share_ar1 = mean(share_ar1,2);
reve_ar1  = mean(reve_ar1,2);
exitrate  = mean(exitrate,2);
entryrate = mean(entryrate,2);
avp       = mean(avp,2);

MM = cat(1,learning,reve_ar1,exitrate,entryrate,avp); 
