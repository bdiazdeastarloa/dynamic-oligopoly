function [MM,summ] = moms(par,p1,x1,y1,solveflag)

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
w0   = [1 1 1 2 M+1 M+1 M+1 M+1]';           % initial industry structure
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
    jumps = 35;                               % # of jumps allowed for per period
    draws = struct();
    draws.any = rand(T*jumps,S);              % draws: any jump 
    draws.ind = rand(T*jumps,S);              % draws: firm-level jumps
    draws.out = rand(T*jumps,S);              % draws: aggregate jump
    eval(['save ',savedraws,' draws']);
end


%% Start simulations.

outcflag = zeros(S,1);
%p_pc5    = zeros(S,1);
%p_pc50   = zeros(S,1);
p_pc90   = zeros(S,1);
p_rat1   = zeros(S,1);
%p_rat2   = zeros(S,1);
enrate   = zeros(S,1);
exrate   = zeros(S,1);
p_ac     = zeros(S,1);
rdar1    = zeros(S,1);
revar1   = zeros(S,1);
%qar1   = zeros(S,1);
lbd_q    = zeros(S,1);
lbd_r    = zeros(S,1);
rd1      = zeros(S,1);
rd2      = zeros(S,1);
summ     = zeros(S,5);
emptyind = zeros(S,1);

if solveflag~=1
    %p_pc5    = ones(S,1)*100;
    %p_pc50   = ones(S,1)*100;
    p_pc90   = ones(S,1)*100;
    p_rat1   = ones(S,1)*100;
    %p_rat2   = ones(S,1)*100;
    enrate   = ones(S,1)*100;
    exrate   = ones(S,1)*100;
    p_ac     = ones(S,1)*100;
    rdar1    = ones(S,1)*100;
    revar1   = ones(S,1)*100;
    %qar1   = ones(S,1)*100;
    lbd_q    = ones(S,1)*100;
    lbd_r    = ones(S,1)*100;
    rd1      = ones(S,1)*100;
    rd2      = ones(S,1)*100;

else   
    for sim=1:S
        drw1 = draws.any(:,sim);
        drw2 = draws.ind(:,sim);
        drw3 = draws.out(:,sim);
        
        try
            % Simulate industry: arrival times and jumps.
            outc = sim_cont(par,p1,x1,y1,di_0,drw1,drw2,drw3,T);
            outcflag(sim,:) = outc.flag;
            emptyind(sim) = outc.empty;
        catch
            display('There was an error while simluating continuous industry');
        end
            
        clear drw1 drw2 drw3
        
        % Fit to discrete time grid (years).
        % AD: aggregate data,
        % [ff,entry,exit,cs,ps,scrap,ecost,p0,q0,sh0,d,w,avprice,c2]
        % PD: panel data.
        [agg,panel] = sim_disc(par,outc,T);
    
        % Define variables for moments
        % Simulated data are sorted by id then period
        nf     = outc.nf;
        id     = panel(:,1);
        const  = ones(size(panel,1),1);

        % ---------------------------------------------------------------------
        % Compute simulated moments.
        % ---------------------------------------------------------------------

        % ---------------------------------------------------------------------
        % Price moments.
        % ---------------------------------------------------------------------
        % Distribution of prices.
        pj  = panel(:,4);
        %p_pc5(sim,:)  = prctile(pj,5);
        %p_pc50(sim,:)  = median(pj);
        p_pc90(sim,:) = prctile(pj,90);
        p_rat1(sim,:) = prctile(pj,90)/prctile(pj,5);
        %p_rat2(sim,:) = prctile(pj,50)/prctile(pj,5);
        
        % Price one-period autocorrelation.
        lagpj =panel(:,15);
        ii = (lagpj>0);
        y  = log(pj(ii,:));
        x  = log(lagpj(ii));
        b  = corrcoef(y,x);
        p_ac(sim,:) = b(2,1);
        
        % ---------------------------------------------------------------------
        % Learning curve.
        % ---------------------------------------------------------------------
        cumq  = panel(:,22);
        cumrd = panel(:,23);
        ii  = (cumq>0);
        iii = (cumrd>=0);
        iii = (iii.*ii==1);
        y  = log(pj(iii,:));

        x = [const(iii) log(cumq(iii)) log(cumrd(iii))];
        [b,~,~] = regress(y,x);
        lbd_q(sim,1) = b(2);
        lbd_r(sim,1) = b(3);
        
        % ---------------------------------------------------------------------
        % Exit and entry rates.
        % ---------------------------------------------------------------------
        entry = agg(:,2)./agg(:,1);
        exit  = agg(:,3)./agg(:,1);
        exrate(sim,1) = nanmean(exit);
        enrate(sim,1) = nanmean(entry);

        % ---------------------------------------------------------------------
        % Revenue targets.
        % ---------------------------------------------------------------------
        % Revenues AR1.
        rev    = panel(:,7);
        lagrev = panel(:,17);
        ii = (lagrev>0);
        y  = log(rev(ii,:));
        x  = log(lagrev(ii));
        [b,~,~] = regress(y,x);
        revar1(sim,:) = b(1);
        
        % ---------------------------------------------------------------------
        % R&D targets.
        % ---------------------------------------------------------------------
        % R&D per unit revenue. 
        lagrd  = panel(:,21);
        rd     = panel(:,5);
        share  = panel(:,10);
        % Generate share change rank.
        meanrd = zeros(nf,1);
        dsh    = zeros(nf,1);
        rdrv = rd./rev;
        for i=1:nf
            sh = share(id==i);
            dsh(i) = sh(end) - sh(1);
            meanrd(i) = mean(rdrv(id==i));
        end
        [~,in] = sort(dsh,1,'descend');
        rd1(sim,:) = meanrd(in(1));
        rd2(sim,:) = meanrd(in(2));
        
        % Average (normalized)  rate of change of RD/Revenue.
%         D = [];
%         for i=1:nf 
%             temp=rdrev(id==i); 
%             temp1=temp(1:sum(id==i)-1); 
%             temp2=temp(2:end); 
%             B=(temp2-temp1)./(temp1+temp2)/2; 
%             D=[D;B]; %#ok
%         end;
%         rdchange = nanmean(D);

        % RD AR1.
        ii = (lagrd>0);
        y  = log(1+rd(ii,:));
        x  = log(1+lagrd(ii));
        [b,~,~] = regress(y,x);
        rdar1(sim,:) = b(1);
        
        % -----------------------------------------------------------------------------
        % Summary: total # of firms, avg # of firms, # of entrants, # of exitors
        % -----------------------------------------------------------------------------
        summ(sim,:) = [nf mean(agg(:,1)) sum(agg(:,2)) sum(agg(:,3)) emptyind(sim)];
         
    end
end

if sum(emptyind)>0
    display(['There were ',num2str(sum(emptyind)),' simulations with an empty industry.']);
end

clear draws agg panel 

%% Compute moments.

% Punish if there was a problem simulating.
in = (outcflag==1);

%p_pc5(in)    = 100;
%p_pc50(in)   = 100;
p_pc90(in)   = 100;
p_rat1(in)   = 100;
%p_rat2(in)   = 100;
enrate(in)   = 100;
exrate(in)   = 100;
p_ac(in)     = 100;
rdar1(in)    = 100;
revar1(in)   = 100;
%qar1(in)     = 100;
lbd_q(in)    = 100;
lbd_r(in)    = 100;
rd1(in)      = 100;
rd2(in)      = 100;

%p_pc5   = mean(p_pc5,1);
%p_pc50  = mean(p_pc50,1); 
p_pc90  = mean(p_pc90,1);
p_rat1  = mean(p_rat1,1);
%p_rat2  = mean(p_rat2,1);
exrate  = nanmean(exrate,1);
enrate  = nanmean(enrate,1);
p_ac    = mean(p_ac,1);
rdar1   = mean(rdar1,1);
revar1  = mean(revar1,1);
%qar1    = mean(qar1,1);
lbd_q   = mean(lbd_q,1);
lbd_r   = mean(lbd_r,1);
rd1     = mean(rd1,1);
rd2     = mean(rd2,1);

clear outcflag 

MM = [
%p_pc5;
%p_pc50;
p_pc90;
p_rat1;
%p_rat2;
enrate;
exrate;
p_ac;
rdar1;
revar1;
%qar1;
lbd_q;
lbd_r;           
rd1;
rd2];

clear p_pc5 p_pc50 p_pc90 p_rat1 p_rat2 enrate exrate p_ac rdar1 revar1 qar1 lbd_q lbd_r rd1 rd2
