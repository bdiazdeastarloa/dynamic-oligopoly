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
p0   = 1;                                    % initial foreign price state
w0   = [1 2 3 3 3 M+1 M+1 M+1]';             % initial industry structure
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
    jumps = 20;                               % # of jumps allowed for per period
    draws = struct();
    draws.any = rand(T*jumps,S);              % draws: any jump 
    draws.ind = rand(T*jumps,S);              % draws: firm-level jumps
    draws.out = rand(T*jumps,S);              % draws: aggregate jump
    eval(['save ',savedraws,' draws']);
end


%% Start simulations.

outcflag = zeros(S,1);
pricepct = zeros(S,1);
revear1  = zeros(S,1);
exrate   = zeros(S,1);
enrate   = zeros(S,1);
rdrev1   = zeros(S,1);
rdrev2   = zeros(S,1);
dp1      = zeros(S,1);
dp2      = zeros(S,1);
drawany  = draws.any;
drawind  = draws.ind;
drawout  = draws.out;

if solveflag~=1
    pricepct = ones(S,1)*100;
    revear1  = ones(S,1)*100;
    exrate   = ones(S,1)*100;
    enrate   = ones(S,1)*100;
    rdrev1   = ones(S,1)*100;
    rdrev2   = ones(S,1)*100;
    dp1      = ones(S,1)*100;
    dp2      = ones(S,1)*100;  
else   
    for sim=1:S
        drw1 = drawany(:,sim);
        drw2 = drawind(:,sim);
        drw3 = drawout(:,sim);
        % Simulate industry: arrival times and jumps.
        outc = sim_cont(par,p1,x1,y1,di_0,drw1,drw2,drw3,T);
        outcflag(sim,:) = outc.flag;
        
        % Fit to discrete time grid (years).
        % AD: aggregate data,
        % [ff,entry,exit,cs,ps,scrap,ecost,p0,q0,sh0,d,w,avprice,c2]
        % PD: panel data.
        [agg,panel] = sim_disc(par,outc,T);

        % Define variables for moments
        % Simulated data are sorted by id then period
        nf     = outc.nf;
        id     = panel(:,1);
        price  = panel(:,4);
        rd     = panel(:,5);
        rev    = panel(:,7);
        share  = panel(:,10);
        lagrev = panel(:,18);

        entry = agg(:,2);
        exit  = agg(:,3);

        const = ones(size(panel,1),1);

        % ---------------------------------------------------------------------
        % Compute simulated moments.
        % ---------------------------------------------------------------------

        % ---------------------------------------------------------------------
        % Leaning curve.
        % ---------------------------------------------------------------------
        % ii = (cumq1>0);
        % ii = (cumq2>0);
        % y  = log(price(ii,:));

        % Generate firm fixed effects matrix
        % temp   = eye(nf);
        % fix_id = zeros(size(y,1),nf); 
        % tempid = id(ii);
        % for i=1:size(y,1)
        %     fix_id(i,:) = (temp(tempid(i),:));
        % end
        % Check for columns of zeros
        % check  = (sum(fix_id)>0);
        % fix_id = fix_id(:,check); 
        % x = [const(ii) log(cumq1(ii)) fix_id(:,2:end)];
        % x = [const(ii) log(cumq2(ii)) fix_id(:,2:end)];

        % [b,~,~] = regress(y,x);
        % rMSE = sqrt(sum(r.^2)/size(y,1));
        % learning(sim,:) = b(2);

        % ---------------------------------------------------------------------
        % Exit and entry rates.
        % ---------------------------------------------------------------------
        exrate(sim,:) = mean(exit);
        enrate(sim,:) = mean(entry);

        % ---------------------------------------------------------------------
        % Log revenues AR1.
        % ---------------------------------------------------------------------
        ii = (lagrev>0);
        y  = log(rev(ii,:));
        x  = [const(ii) log(lagrev(ii))];
        [b,~,~] = regress(y,x);
        % rMSE    = sqrt(sum(r.^2)/(size(y,1)-2));
        revear1(sim,:) = b(2);

        % ---------------------------------------------------------------------
        % Percentile of the price distribution.
        % ---------------------------------------------------------------------
        pricepct(sim,:) = prctile(price,10);

        % ---------------------------------------------------------------------
        % Mean R&D per unit revenue.
        % ---------------------------------------------------------------------
        % Generate share change rank.
        meanrd = zeros(nf,1);
        dsh    = zeros(nf,1);
        dp     = zeros(nf,1);
        temp = rd./rev;
        for i=1:nf
            sh = share(id==i);
            p  = price(id==i);
            dsh(i) = sh(end) - sh(1);
            dp(i)  = log(p(end)/p(1));
            meanrd(i) = mean(temp(id==i));
        end
        [~,in] = sort(dsh,1,'descend');
        rdrev1(sim,:) = meanrd(in(1));
        rdrev2(sim,:) = meanrd(in(2));
        dp1(sim,:)    = dp(in(1));
        dp2(sim,:)    = dp(in(2));
    end
end

%% Compute moments.

% Punish if there was a problem simulating.
in = (outcflag==1);

pricepct(in) = 100;
revear1(in) = 100;
exrate(in) = 100;
enrate(in) = 100;
rdrev1(in) = 100;
rdrev2(in) = 100;
dp1(in) = 100;
dp2(in) = 100;

pricepct = mean(pricepct,1);
revear1  = mean(revear1,1);
exrate   = mean(exrate,1);
enrate   = mean(enrate,1);
rdrev1   = mean(rdrev1,1);
rdrev2   = mean(rdrev2,1);
dp1      = mean(dp1,1);
dp2      = mean(dp2,1);

MM = cat(1,pricepct,revear1,exrate,enrate,rdrev1,rdrev2,dp1,dp2); 
