function MM = moms(par,p1,x1,y1,solveflag)

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
%
% moms:
% - creates moments to match data.
% - sim_cont: simulates arrival times and outcomes.
% - sim_disc: fits simulation in a discrete time grid. 
%
% Written by Bernardo Diaz de Astarloa @ PSU 2014.
% =========================================================================

%% Simulation settings and parameters.

p0size = par.p0size;
M      = par.M;

mkt0 = 1;                                    % initial demand state
p0   = p0size;                               % initial foreign price state
w0   = [1 1 1 2 4 M+1 M+1 M+1]';             % initial industry structure
d0   = p0+(mkt0-1)*p0size;                   % index of exogenous shock
di_0 = [d0;w0];

T = 11;                                      % # of periods to simulate
S = 500;                                     % # of simulations to run

% Random draws for simulations (#jumps,#simulations)
savedraws = sprintf('draws_T%dS%d.mat',[T,S]); 
if exist(savedraws,'file')
    load(savedraws)
else
    rng(8282);
    jumps = 30;                               % # of jumps allowed for per period
    draws = struct();
    draws.any = rand(T*jumps,S);              % draws: any jump 
    draws.ind = rand(T*jumps,S);              % draws: firm-level jumps
    draws.out = rand(T*jumps,S);              % draws: aggregate jump
    eval(['save ',savedraws,' draws']);
end


%% Start simulations.

learning  = zeros(S,1);
share_ar1 = zeros(S,1);
reve_ar1  = zeros(S,1);
exitrate  = zeros(S,1);
entryrate = zeros(S,1);
avp       = zeros(S,1);

if solveflag~=1
    learning  = ones(S,1)*100;
    share_ar1 = ones(S,1)*100;
    reve_ar1  = ones(S,1)*100;
    exitrate  = ones(S,1)*100;
    entryrate = ones(S,1)*100;
    avp       = ones(S,1)*100;
else   
    for sim=1:S
        % Simulate industry: arrival times and jumps.
        outc = sim_cont(par,p1,x1,y1,di_0,draws,T,sim);
        
        if outc.flag==1
            disp('Something went wrong in simulation.');
            learning  = ones(S,1)*100;
            share_ar1 = ones(S,1)*100;
            reve_ar1  = ones(S,1)*100;
            exitrate  = ones(S,1)*100;
            entryrate = ones(S,1)*100;
            avp       = ones(S,1)*100;
            break
        else   
            % Fit to discrete time grid (years).
            % AD: aggregate data | PD: panel data.
            [agg,panel,~] = sim_disc(par,outc,T);

            % Define variables for moments
            % Simulated data are sorted by id then period
            nf     = outc.nf;
            id     = panel(:,1);
            year   = panel(:,2);
            price  = panel(:,4);
            rev    = panel(:,6);
            share  = panel(:,9);

            lagrev   = panel(:,17);
            lagshare = panel(:,18);
            
            cumq1    = panel(:,14);
            cumq2    = panel(:,21);
            
            entry = agg(:,2);
            exit  = agg(:,3);
            
            const = ones(size(panel,1),1);

            % ---------------------------------------------------------------------
            % Compute simulated moments
            % ---------------------------------------------------------------------
            
            % ---------------------------------------------------------------------
            % Leaning curve
            % ---------------------------------------------------------------------
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
            learning(sim,:) = b(2);

            
            % ---------------------------------------------------------------------
            % Exit and entry rates
            % ---------------------------------------------------------------------
            exitrate(sim,:) = mean(exit);
            entryrate(sim,:) = mean(entry);
            
            
            % ---------------------------------------------------------------------
            % Log shares AR1
            % ---------------------------------------------------------------------
            ii = (lagshare>0);
            y  = log(share(ii,:));
            x  = [const(ii) year(ii) log(lagshare(ii))];
            [b,~,~] = regress(y,x);
            % rMSE    = sqrt(sum(r.^2)/(size(y,1)-2));
            share_ar1(sim,:) = b(3);
            
            
            % ---------------------------------------------------------------------
            % Log revenues AR1
            % ---------------------------------------------------------------------
            ii = (lagrev>0);
            y  = log(rev(ii,:));
            x  = [const(ii) year(ii) log(lagrev(ii))];
            [b,~,~] = regress(y,x);
            % rMSE    = sqrt(sum(r.^2)/(size(y,1)-2));
            reve_ar1(sim,:) = b(3);

            
            % ---------------------------------------------------------------------
            % Mean price
            % ---------------------------------------------------------------------
            avp(sim,:) = mean(price(price>0));
            
        end
    end
end


%% Compute moments.

learning  = mean(learning,1);
share_ar1 = mean(share_ar1,1);
reve_ar1  = mean(reve_ar1,1);
exitrate  = mean(exitrate,1);
entryrate = mean(entryrate,1);
avp       = mean(avp,1);

MM = cat(1,learning,reve_ar1,exitrate,entryrate,avp); 
