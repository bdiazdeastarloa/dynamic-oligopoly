function [m,sd] = sim_outcomes(P,X,E,par)

%% Industry Simulations
T = par.T;              % # of periods to simulate
S = par.S;              % # of simulations to run

wstart  = par.wstart;   % initial industry structure
w0start = par.w0start;  % starting outside good price

nfirms   = zeros(T,S);
entrants = zeros(T,S);
exitors  = zeros(T,S);
csurplus = zeros(T,S);
profits  = zeros(T,S);
scrap    = zeros(T,S);
ecost    = zeros(T,S);
p0       = zeros(T,S);
q0       = zeros(T,S);
share0   = zeros(T,S);
avprice  = zeros(T,S);
c2       = zeros(T,S);
omega    = zeros(T,S);

parfor s=1:S  
    % Simulate industry: arrival times and jumps
    [nf,times_c,~,exit_c,entry_c,fprod_c,pij_c,~,p0_c,q0_c,mgn_c,quant_c,rev_c,ecostj_c,scrapj_c,ecost_c,scrap_c,csurplus_c,aggprof_c,nfirms_c,simflag] = sim_cont(par,T,P,E,X,wstart,w0start,s);

    if simflag==1
        disp('Something went wrong in simulation');
        nfirms(:,s)   = NaN(T,1);
        entrants(:,s) = NaN(T,1);
        exitors(:,s)  = NaN(T,1);
        csurplus(:,s) = NaN(T,1);
        profits(:,s)  = NaN(T,1);
        scrap(:,s)    = NaN(T,1);
        ecost(:,s)    = NaN(T,1);
        p0(:,s)       = NaN(T,1);
        q0(:,s)       = NaN(T,1);
        share0(:,s)   = NaN(T,1);
        avprice(:,s)  = NaN(T,1);
        c2(:,s)       = NaN(T,1);
        omega(:,s)    = NaN(T,1);
    else   
        % Fit to discrete time grid (years)
        % Aggregate data & Panel data
        % Simulated panel data are sorted by id then period
        [agg,~,avg] = sim_disc(par,nf,times_c,exit_c,entry_c,fprod_c,pij_c,p0_c,q0_c,mgn_c,quant_c,rev_c,ecostj_c,scrapj_c,ecost_c,scrap_c,csurplus_c,aggprof_c,nfirms_c,T);

        % -----------------------------------------------------------------
        % Aggregate data
        nfirms(:,s)   = agg(:,1);
        entrants(:,s) = agg(:,2);
        exitors(:,s)  = agg(:,3);
        csurplus(:,s) = agg(:,4);
        profits(:,s)  = agg(:,5);
        scrap(:,s)    = agg(:,6);
        ecost(:,s)    = agg(:,7);
        p0(:,s)       = agg(:,8);
        q0(:,s)       = agg(:,9);
        share0(:,s)   = agg(:,10);

        % -----------------------------------------------------------------
        % Averaged data
        omega(:,s)   = avg(:,1);
        avprice(:,s) = avg(:,2);
        c2(:,s)      = avg(:,3);
    end
end

m.nfirms   = nanmean(nfirms,2);
m.entrants = nanmean(entrants,2);
m.exitors  = nanmean(exitors,2);
m.csurplus = nanmean(csurplus,2);
m.profits  = nanmean(profits,2);
m.scrap    = nanmean(scrap,2);
m.ecost    = nanmean(ecost,2);
m.p0       = nanmean(p0,2);
m.q0       = nanmean(q0,2);
m.share0   = nanmean(share0,2);
m.avprice  = nanmean(avprice,2);
m.c2       = nanmean(c2,2);
m.omega    = nanmean(omega,2);

sd.nfirms   = nanstd(nfirms,0,2);
sd.entrants = nanstd(entrants,0,2);
sd.exitors  = nanstd(exitors,0,2);
sd.csurplus = nanstd(csurplus,0,2);
sd.profits  = nanstd(profits,0,2);
sd.scrap    = nanstd(scrap,0,2);
sd.ecost    = nanstd(ecost,0,2);
sd.p0       = nanstd(p0,0,2);
sd.q0       = nanstd(q0,0,2);
sd.share0   = nanstd(share0,0,2);
sd.avprice  = nanstd(avprice,0,2);
sd.concentr = nanstd(c2,0,2);
sd.omega    = nanstd(omega,0,2);


