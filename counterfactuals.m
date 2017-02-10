function counterfactuals

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
% 
% Simulate industry paths under alternative scenarios and calculate
% indicators.
%
% Output is saved to MAT-files (one for each scenario) that can be used in 
% 'count_out.m'.
%
% Written by Bernardo Diaz de Astarloa @ PSU March 2015
% =========================================================================

clear all
format long 

% Compute equilibria. 
%{
% Calibrated parameter value.
parguess = [0.85038     0.457336     0.267624   0.00390264        6.121      5.15614     0.367204      11.7505     0.709154];

for countfact=0:1
    
    alpha1  = parguess(1);              % r&d hazard: scale
    alpha2  = parguess(2);              % r&d hazard: elasticity
    eta1    = parguess(3);              % learning hazard: scale
    eta2    = parguess(4);              % learning hazard: elasticity
    etax1   = parguess(5);              % exit opportunity hazard
    etax2   = 0;                        % exit opportunity hazard
    etae1   = parguess(6);              % entry opportunity hazard
    etae2   = 0 ;                       % entry opportunity hazard
    delta   = parguess(7);              % negative productivity shock hazard
    c0      = parguess(8);              % cost function level
    beta    = parguess(9);              % cost function curvature
    Phi_hi  = 100;                      % scrap value upper bound
    Phi_lo  = 50;                       % scrap value lower bound
    Phie_hi = 150;                      % entry cost upper bound
    Phie_lo = 100;                      % entry cost lower bound
    alpha   = 2.779;                    % price coefficient
    
    % Set counterfactual scenario.
    counter = countfact; 
    % Assign parameters.
    setpar;
    counteq = ['results/count_c' int2str(counter) '.mat'];
    
    % Compute initial values for policies and value function.
    initfile2 = ['initials_a' sprintf('%5.4f', par.alpha) '.mat'];
    sol_last  = 'lastloop_sol.mat';
    
    if exist(sol_last,'file')
        % Load converged values from last parameter guess.
        load(sol_last)
        V0 = V1;                            
        p0 = p1;                            
        x0 = x1;                            
        y0 = y1;                            
        clear V1 p1 x1 y1 info_
    elseif exist(initfile2, 'file')
        % Load static equilibrium if it exists.
        load(initfile2)
    else
        % Set initial prices and values to static game solution
        [p0,V0] = init(par);                
        x0 = zeros(N,S);                    
        y0 = zeros(N,S);                    
        eval(['save ',initfile2,' p0 V0 x0 y0']);
    end
    % Assign initial values to par structure.
    par.V0 = V0;
    par.x0 = x0;
    par.p0 = p0;
    par.y0 = y0;
    clear V0 x0 p0 y0

    % Solve for a MPE.
    [info_,V1,x1,p1,y1] = compMPE(par);

    solveflag = info_(1);
    if solveflag==1
        existence(par,V1,p1);
        save(counteq,'V1','p1','x1','y1','par');
        save(sol_last,'V1','p1','x1','y1','info_');        
    else
        disp('WARNING: something went wrong when solving counterfactual.')
    end
end
%}


%% Simulate counterfactuals.
for counter = 0:1;
    [simcount_m,simcount_sd] = sim_counter(counter);          %#ok
    simcount = ['results/simcount_c' int2str(counter) '.mat'];
    save(simcount,'simcount_m','simcount_sd');
end

% Counterfactuals output.
count_out;

end


%% Simulate counterfactuals.
function [m,sd] = sim_counter(counter)

TT = 25;                                     % # of periods to simulate
T = 15;                                      % # of periods to simulate before shock
S = 500;                                     % # of simulations to run

% Load baseline.
load(['results/count_c' int2str(0) '.mat']);
parpre = par;
p1pre  = p1;
x1pre  = x1;
y1pre  = y1;

M   = par.M;
N   = par.N;
D   = par.D;
rho = par.rho;
p0size = par.p0size;

% Initial conditions.
mkt_in = 1;                                    % initial demand state
p0_in  = 1;                                    % initial foreign price state
w_in   = [1 1 1 2 M+1 M+1 M+1 M+1]';           % initial industry structure
d_in   = p0_in+(mkt_in-1)*p0size;              % index of exogenous shock
di_in  = [d_in;w_in];

clear par p1 x1 y1

% Load post-baseline (indexed by 'counter').
load(['results/count_c' int2str(counter) '.mat']);
parpost = par;
p1post  = p1;
x1post  = x1;
y1post  = y1;
clear par p1 x1 y1

% Random draws for simulations (#jumps,#simulations)
savedrawspre = sprintf('drawspre_T%dS%d.mat',[T,S]); 
if exist(savedrawspre,'file')
    load(savedrawspre)
else
    rng(2282);
    jumps = 40;                                      % # of jumps allowed for per period
    drawspre = struct();
    drawspre.any = rand(T*jumps,S);              % draws: any jump 
    drawspre.ind = rand(T*jumps,S);              % draws: firm-level jumps
    drawspre.out = rand(T*jumps,S);              % draws: aggregate jump
    eval(['save ',savedrawspre,' drawspre']);
end
load(savedrawspre);
savedrawspost = sprintf('drawspost_T%dS%d.mat',[TT-T,S]); 
if exist(savedrawspost,'file')
    load(savedrawspost)
else
    rng(2212);
    jumps = 40;                                      % # of jumps allowed for per period
    drawspost = struct();
    drawspost.any = rand((TT-T)*jumps,S);              % draws: any jump 
    drawspost.ind = rand((TT-T)*jumps,S);              % draws: firm-level jumps
    drawspost.out = rand((TT-T)*jumps,S);              % draws: aggregate jump
    eval(['save ',savedrawspost,' drawspost']);
end

nfirms   = zeros(TT,S);
entry    = zeros(TT,S);
exit     = zeros(TT,S);
csurplus = zeros(TT,S);
profits  = zeros(TT,S);
netprof  = zeros(TT,S);
scrap    = zeros(TT,S);
ecost    = zeros(TT,S);
pfor     = zeros(TT,S);
share0   = zeros(TT,S);
omega    = zeros(TT,S);
avprice  = zeros(TT,S);
c2       = zeros(TT,S);
avrd     = zeros(TT,S);
avmgc    = zeros(TT,S);
rdjump   = zeros(TT,S);
lbdjump  = zeros(TT,S);
d_cs     = zeros(S,1);
d_psx    = zeros(S,1);
d_sc     = zeros(S,1);
d_ec     = zeros(S,1);

for s=1:S
    % display(['Simulation ',int2str(s)]);
    % Simulate industry pre-shock.
    drw1pre = drawspre.any(:,s);
    drw2pre = drawspre.ind(:,s);
    drw3pre = drawspre.out(:,s);
    outc = sim_cont(parpre,p1pre,x1pre,y1pre,di_in,drw1pre,drw2pre,drw3pre,T);
    
    if outc.flag==1
        disp('Something went wrong in simulation (pre-shock).');
        nfirms(1:TT,s)   = NaN(TT,1);
        entry(1:TT,s)    = NaN(TT,1);
        exit(1:TT,s)     = NaN(TT,1);
        csurplus(1:TT,s) = NaN(TT,1);
        profits(1:TT,s)  = NaN(TT,1);
        netprof(1:TT,s)  = NaN(TT,1);
        scrap(1:TT,s)    = NaN(TT,1);
        ecost(1:TT,s)    = NaN(TT,1);
        pfor(1:TT,s)     = NaN(TT,1);
        share0(1:TT,s)   = NaN(TT,1);
        omega(1:TT,s)    = NaN(TT,1);
        avprice(1:TT,s)  = NaN(TT,1);
        c2(1:TT,s)       = NaN(TT,1);
        avrd(1:TT,s)     = NaN(TT,1);
        avmgc(1:TT,s)    = NaN(TT,1);
        rdjump(1:TT,s)   = NaN(TT,1);
        lbdjump(1:TT,s)  = NaN(TT,1);
        d_cs(s)          = NaN(1,1);
        d_psx(s)    	 = NaN(1,1);
        d_sc(s)          = NaN(1,1);
        d_ec(s)          = NaN(1,1);
        
    else   
        % Fit to discrete time grid (years).
        [agg,~] = sim_disc(parpre,outc,T);
        
        nfirms(1:T,s)   = agg(:,1);
        entry(1:T,s)    = agg(:,2);
        exit(1:T,s)     = agg(:,3);
        csurplus(1:T,s) = agg(:,4);
        profits(1:T,s)  = agg(:,5);
        netprof(1:T,s)  = agg(:,6);
        scrap(1:T,s)    = agg(:,7);
        ecost(1:T,s)    = agg(:,8);
        pfor(1:T,s)     = agg(:,9);
        share0(1:T,s)   = agg(:,11);
        omega(1:T,s)    = agg(:,13);
        avprice(1:T,s)  = agg(:,14);
        c2(1:T,s)       = agg(:,15);
        avrd(1:T,s)     = agg(:,16);
        avmgc(1:T,s)    = agg(:,18);
        rdjump(1:T,s)   = agg(:,19);
        lbdjump(1:T,s)  = agg(:,20);
        
        % Collect welfare objects for later.
        welf_pre = [outc.times outc.cs outc.psx outc.scrap outc.ecost];
        
        % Simulate industry post-shock.
        drw1post = drawspost.any(:,s);
        drw2post = drawspost.ind(:,s);
        drw3post = drawspost.out(:,s);

        w0 = outc.fprod(end,(outc.fprod(end,:)>0));
        w0 = sort(w0); 
        w0 = [w0 (M+1)*ones(1,N-size(w0,2))]';
        mkt0 = 1;                                    % initial demand state
        p0   = find(parpre.pfor==outc.p0(end));      % initial foreign price state
        d0   = p0+(mkt0-1)*p0size;                   % index of exogenous shock
        di_0 = [d0;w0];

        outc = sim_cont(parpost,p1post,x1post,y1post,di_0,drw1post,drw2post,drw3post,TT-T);
        if outc.flag==1
            disp('Something went wrong in simulation (post-shock).');
            nfirms(1:TT,s)   = NaN(TT,1);
            entry(1:TT,s)    = NaN(TT,1);
            exit(1:TT,s)     = NaN(TT,1);
            csurplus(1:TT,s) = NaN(TT,1);
            profits(1:TT,s)  = NaN(TT,1);
            netprof(1:TT,s)  = NaN(TT,1);
            scrap(1:TT,s)    = NaN(TT,1);
            ecost(1:TT,s)    = NaN(TT,1);
            pfor(1:TT,s)     = NaN(TT,1);
            share0(1:TT,s)   = NaN(TT,1);
            omega(1:TT,s)    = NaN(TT,1);
            avprice(1:TT,s)  = NaN(TT,1);
            c2(1:TT,s)       = NaN(TT,1);
            avrd(1:TT,s)     = NaN(TT,1);
            avmgc(1:TT,s)    = NaN(TT,1);
            rdjump(1:TT,s)   = NaN(TT,1);
            lbdjump(1:TT,s)  = NaN(TT,1);
            d_cs(s)          = NaN(1,1);
            d_psx(s)    	 = NaN(1,1);
            d_sc(s)          = NaN(1,1);
            d_ec(s)          = NaN(1,1);
        else
            % Fit to discrete time grid (years).
            [agg,~] = sim_disc(parpost,outc,TT-T);

            % Aggregate data
            nfirms(T+1:end,s)   = agg(:,1);
            entry(T+1:end,s)    = agg(:,2);
            exit(T+1:end,s)     = agg(:,3);
            csurplus(T+1:end,s) = agg(:,4);
            profits(T+1:end,s)  = agg(:,5);
            netprof(T+1:end,s)  = agg(:,6);
            scrap(T+1:end,s)    = agg(:,7);
            ecost(T+1:end,s)    = agg(:,8);
            pfor(T+1:end,s)     = agg(:,9);
            share0(T+1:end,s)   = agg(:,11);
            omega(T+1:end,s)    = agg(:,13);
            avprice(T+1:end,s)  = agg(:,14);
            c2(T+1:end,s)       = agg(:,15);
            avrd(T+1:end,s)     = agg(:,16);
            avmgc(T+1:end,s)    = agg(:,18);
            rdjump(T+1:end,s)   = agg(:,19);
            lbdjump(T+1:end,s)  = agg(:,20);

            % Collect welfare objects for later.
            welf_post = [outc.times+T outc.cs outc.psx outc.scrap outc.ecost];
            
            % Compute welfare measures.
            [d_cs(s),d_psx(s),d_sc(s),d_ec(s)] =  welfare(rho,welf_pre,welf_post,TT);
        end
    end  
end

m.nfirms   = nanmean(nfirms,2);
m.entry    = nanmean(entry,2);
m.exit     = nanmean(exit,2);
m.csurplus = nanmean(csurplus,2);
m.profits  = nanmean(profits,2);
m.netprofits  = nanmean(netprof,2);
m.scrap    = nanmean(scrap,2);
m.ecost    = nanmean(ecost,2);
m.pfor     = nanmean(pfor,2);
m.share0   = nanmean(share0,2);
m.omega    = nanmean(omega,2);
m.avprice  = nanmean(avprice,2);
m.c2       = nanmean(c2,2);
m.avrd     = nanmean(avrd,2);
m.avmgc    = nanmean(avmgc,2);
m.rdjump   = nanmean(rdjump,2);
m.lbdjump  = nanmean(lbdjump,2);
m.d_cs     = nanmean(d_cs,1);
m.d_psx    = nanmean(d_psx,1);
m.d_sc     = nanmean(d_sc,1);
m.d_ec     = nanmean(d_ec,1);

sd.nfirms   = nanstd(nfirms,0,2);
sd.entry    = nanstd(entry,0,2);
sd.exit     = nanstd(exit,0,2);
sd.csurplus = nanstd(csurplus,0,2);
sd.profits  = nanstd(profits,0,2);
sd.netprofits  = nanstd(netprof,0,2);
sd.scrap    = nanstd(scrap,0,2);
sd.ecost    = nanstd(ecost,0,2);
sd.pfor     = nanstd(pfor,0,2);
sd.share0   = nanstd(share0,0,2);
sd.omega    = nanstd(omega,0,2);
sd.avprice  = nanstd(avprice,0,2);
sd.c2       = nanstd(c2,0,2);
sd.avrd     = nanstd(avrd,0,2);
sd.avmgc    = nanstd(avmgc,0,2);
sd.rdjump   = nanstd(rdjump,0,2);
sd.lbdjump  = nanstd(lbdjump,0,2);
sd.d_cs     = nanstd(d_cs,0,1);
sd.d_psx    = nanstd(d_psx,0,1);
sd.d_sc     = nanstd(d_sc,0,1);
sd.d_ec     = nanstd(d_ec,0,1);

end


%% Compute welfare measure.
function [cs,psx,scrap,entry] = welfare(rho,welf_pre,welf_post,TT)

time  = [welf_pre(:,1);welf_post(:,1)];
cs    = [welf_pre(:,2);welf_post(:,2)];
psx   = [welf_pre(:,3);welf_post(:,3)];
scrap = [welf_pre(:,4);welf_post(:,4)];
entry = [welf_pre(:,5);welf_post(:,5)];

disc1 = exp(-rho*time);
disc2 = [exp(-rho*time(2:end)) ; exp(-rho*TT)];
disc  = (disc2 - disc1)/(-rho);

cs    = sum(disc.*cs);
psx   = sum(disc.*psx);
scrap = sum(disc1.*scrap);
entry = sum(disc1.*entry);

end
