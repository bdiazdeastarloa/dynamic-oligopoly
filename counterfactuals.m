function counterfactuals

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

clear all
format long

% Compute equilibria. 
%{
for countfact=0:1
    % Calibrated parameter value.
    parguess = [0.131   1.0900   2.500   2.500   11.1   0.7];

    alpha2  = parguess(1);              %#ok r&d hazard: elasticity 
    alpha1  = 1;                        %#ok r&d hazard: scale
    eta1    = parguess(2);              %#ok learning hazard: scale
    eta2    = 0;                        %#ok learning hazard: elasticity
    etax1   = parguess(3);              %#ok exit opportunity hazard
    etax2   = 0;                        %#ok exit opportunity hazard
    etae1   = parguess(4);              %#ok entry opportunity hazard
    etae2   = 0 ;                       %#ok entry opportunity hazard
    c0      = parguess(5);              %#ok cost function level
    delta   = parguess(6);              %#ok negative productivity shock hazard
    Phi_hi  = 100;                      %#ok scrap value upper bound
    Phi_lo  = 50;                       %#ok scrap value lower bound
    Phie_hi = 150;                      %#ok entry cost upper bound
    Phie_lo = 100;                      %#ok entry cost lower bound
    alpha   = 2.779;                    % price coefficient
    
    % Set counterfactual scenario.
    counter = countfact; 
    % Assign parameters.
    setpar;
    counteq = ['results/count_N' int2str(N) 'M' int2str(M) 'D' int2str(D) '_c' int2str(counter) '.mat'];
    
    % Compute initial values for policies and value function.
    initfile2 = ['initials_N' int2str(N) 'M' int2str(M) 'D' int2str(D) 'a' sprintf('%5.4f', alpha) '.mat'];
    sol_last  = ['lastloop_N' int2str(N) 'M' int2str(M) 'D' int2str(D) '.mat'];
    
    if exist(sol_last,'file')
        % Load converged values from last parameter guess.
        load(sol_last)
        V0 = V1;                            %#ok
        p0 = p1;                            %#ok
        x0 = x1;                            %#ok
        y0 = y1;                            %#ok
        clear V1 p1 x1 y1
    elseif exist(initfile2, 'file')
        % Load static equilibrium if it exists.
        load(initfile2)
    else
        % Set initial prices and values to static game solution
        [p0,V0] = init(par);                %#ok
        x0 = zeros(N,S);                    %#ok
        y0 = zeros(N,S);                    %#ok
        eval(['save ',initfile2,' p0 V0 x0 y0']);
    end
    % Save variables and clear them from workspace.
    eval(['save ',mexinp,vars]);
    eval(['clear ',vars]);

    % Solve for a MPE.
    compMPE(mexinp);
    load(mexinp);

    solveflag = info_(1);
    if solveflag==1
        save(counteq,'V1','p1','x1','y1','par');
        save(sol_last,'V1','p1','x1','y1','info_');        
    else
        disp('WARNING: something went wrong when solving counterfactual.')
    end
    
    clearvars -except N M D 
end
%}

N = 8;
M = 6;
D = 5;

% Simulate counterfactuals.
for counter = 0:1;
    [simcount_m,simcount_sd] = sim_counter(N,M,D,counter);          %#ok
    simcount = ['results/simcount_N' int2str(N) 'M' int2str(M) 'D' int2str(D) '_c' int2str(counter) '.mat'];
    save(simcount,'simcount_m','simcount_sd','agg');
end

% Figures.
count_fig;

end


%% Simulate counterfactuals.
function [m,sd] = sim_counter(N,M,D,counter)

% Load baseline.
load(['results/count_N' int2str(N) 'M' int2str(M) 'D' int2str(D) '_c' int2str(0) '.mat']);
parpre = par;
p1pre  = p1;
x1pre  = x1;
y1pre  = y1;

p0size = par.p0size;
M   = par.M;
N   = par.N;
D   = par.D;
rho = par.rho;
mkt0 = 1;                                    % initial demand state
p0   = 1;                                    % initial foreign price state
w0   = [1 2 3 3 3 M+1 M+1 M+1]';             % initial industry structure
d0   = p0+(mkt0-1)*p0size;                   % index of exogenous shock
di_0 = [d0;w0];

TT = 20;                                     % # of periods to simulate
T = 11;                                      % # of periods to simulate before shock
S = 500;                                     % # of simulations to run

clear par p1 x1 y1
load(['results/count_N' int2str(N) 'M' int2str(M) 'D' int2str(D) '_c' int2str(counter) '.mat']);
parpost = par;
p1post  = p1;
x1post  = x1;
y1post  = y1;
clear par p1 x1 y1

% Random draws for simulations (#jumps,#simulations)
savedraws = sprintf('draws_T%dS%d.mat',[T,S]); 
load(savedraws);
savedrawspost = sprintf('drawspost_T%dS%d.mat',[TT-T,S]); 
if exist(savedrawspost,'file')
    load(savedrawspost)
else
    rng(2212);
    jumps = 20;                                      % # of jumps allowed for per period
    drawspost = struct();
    drawspost.any = rand((TT-T)*jumps,S);              % draws: any jump 
    drawspost.ind = rand((TT-T)*jumps,S);              % draws: firm-level jumps
    drawspost.out = rand((TT-T)*jumps,S);              % draws: aggregate jump
    eval(['save ',savedrawspost,' drawspost']);
end

nfirms   = zeros(TT,S);
csurplus = zeros(TT,S);
profits  = zeros(TT,S);
scrap    = zeros(TT,S);
ecost    = zeros(TT,S);
pfor     = zeros(TT,S);
share0   = zeros(TT,S);
avprice  = zeros(TT,S);
c2       = zeros(TT,S);
avrd     = zeros(TT,S);
avmgc    = zeros(TT,S);

for s=1:S      
    % Simulate industry pre-shock.
    drw1pre = draws.any(:,s);
    drw2pre = draws.ind(:,s);
    drw3pre = draws.out(:,s);
    outc = sim_cont(parpre,p1pre,x1pre,y1pre,di_0,drw1pre,drw2pre,drw3pre,T);

    if outc.flag==1
        disp('Something went wrong in simulation');
        nfirms(1:T,s)   = NaN(T,1);
        csurplus(1:T,s) = NaN(T,1);
        profits(1:T,s)  = NaN(T,1);
        scrap(1:T,s)    = NaN(T,1);
        ecost(1:T,s)    = NaN(T,1);
        pfor(1:T,s)     = NaN(T,1);
        share0(1:T,s)   = NaN(T,1);
        avprice(1:T,s)  = NaN(T,1);
        c2(1:T,s)       = NaN(T,1);
        avrd(1:T,s)     = NaN(T,1);
        avmgc(1:T,s)    = NaN(T,1);
    else   
        % Fit to discrete time grid (years).
        discount = exp(-rho*(outc.times));
        outc.cs = outc.cs.*discount;
        outc.ps = outc.ps.*discount;
        [agg,~] = sim_disc(parpre,outc,T);
        
        nfirms(1:T,s)   = agg(:,1);
        csurplus(1:T,s) = agg(:,4);
        profits(1:T,s)  = agg(:,5);
        scrap(1:T,s)    = agg(:,6);
        ecost(1:T,s)    = agg(:,7);
        pfor(1:T,s)     = agg(:,8);
        share0(1:T,s)   = agg(:,10);
        avprice(1:T,s)  = agg(:,13);
        c2(1:T,s)       = agg(:,14);
        avrd(1:T,s)     = agg(:,15);
        avmgc(1:T,s)    = agg(:,17);
    end
    
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
        disp('Something went wrong in simulation');
        nfirms(TT-T+1:end,s)   = NaN(TT-T,1);
        csurplus(TT-T+1:end,s) = NaN(TT-T,1);
        profits(TT-T+1:end,s)  = NaN(TT-T,1);
        scrap(TT-T+1:end,s)    = NaN(TT-T,1);
        ecost(TT-T+1:end,s)    = NaN(TT-T,1);
        pfor(TT-T+1:end,s)     = NaN(TT-T,1);
        share0(TT-T+1:end,s)   = NaN(TT-T,1);
        avprice(TT-T+1:end,s)  = NaN(TT-T,1);
        c2(TT-T+1:end,s)       = NaN(TT-T,1);
        avrd(TT-T+1:end,s)     = NaN(TT-T,1);
        avmgc(TT-T+1:end,s)    = NaN(TT-T,1);
    else    
        % Fit to discrete time grid (years).
        discount = exp(-rho*(T+outc.times));
        outc.cs = outc.cs.*discount;
        outc.ps = outc.ps.*discount;
        [agg,~] = sim_disc(parpost,outc,TT-T);
        
        % Aggregate data
        nfirms(T+1:end,s)   = agg(:,1);
        csurplus(T+1:end,s) = agg(:,4);
        profits(T+1:end,s)  = agg(:,5);
        scrap(T+1:end,s)    = agg(:,6);
        ecost(T+1:end,s)    = agg(:,7);
        pfor(T+1:end,s)     = agg(:,8);
        share0(T+1:end,s)   = agg(:,10);
        avprice(T+1:end,s)  = agg(:,13);
        c2(T+1:end,s)       = agg(:,14);
        avrd(T+1:end,s)     = agg(:,15);
        avmgc(T+1:end,s)    = agg(:,17);
    end
end

m.csurplus = nanmean(csurplus,2);
m.profits  = nanmean(profits,2);
m.scrap    = nanmean(scrap,2);
m.ecost    = nanmean(ecost,2);
m.pfor     = nanmean(pfor,2);
m.share0   = nanmean(share0,2);
m.avprice  = nanmean(avprice,2);
m.c2       = nanmean(c2,2);
m.avrd     = nanmean(avrd,2);
m.avmgc    = nanmean(avmgc,2);

sd.csurplus = nanstd(csurplus,0,2);
sd.profits  = nanstd(profits,0,2);
sd.scrap    = nanstd(scrap,0,2);
sd.ecost    = nanstd(ecost,0,2);
sd.pfor     = nanstd(pfor,0,2);
sd.share0   = nanstd(share0,0,2);
sd.avprice  = nanstd(avprice,0,2);
sd.c2       = nanstd(c2,0,2);
sd.avrd     = nanstd(avrd,0,2);
sd.avmgc    = nanstd(avmgc,0,2);

end


%% Couterfatuals figures.
function count_fig

load results/simcount_N8M6D5_c0.mat
avp_0 = simcount_m.avprice;
p0_0      = simcount_m.pfor;
share0_0  = simcount_m.share0;
c2_0      = simcount_m.c2;
csurp_0   = simcount_m.csurplus;
profits_0 = simcount_m.profits;
avrd_0    = simcount_m.avrd;
avmgc_0   = simcount_m.avmgc;

load results/simcount_N8M6D5_c1.mat
avp_1 = simcount_m.avprice;
p0_1      = simcount_m.pfor;
share0_1  = simcount_m.share0;
c2_1      = simcount_m.c2;
csurp_1   = simcount_m.csurplus;
profits_1 = simcount_m.profits;
avrd_1    = simcount_m.avrd;
avmgc_1   = simcount_m.avmgc;

% Domestic average selling price
figure('visible','off')
plot(avp_0,'-','LineWidth',3);
hold on
plot(avp_1,'r--','LineWidth',3);
xlabel('Period','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
ylabel('$/watt','FontSize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_avprice.pdf')


% Imports price
figure('visible','off')
plot(p0_0,'-','LineWidth',3);
hold on
plot(p0_1,'r--','LineWidth',3);
xlabel('Period','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
ylabel('$/Watt','FontSize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_p0.pdf')

% Imports share
figure('visible','off')
plot(share0_0,'-','LineWidth',2);
hold on
plot(share0_1,'r--','LineWidth',2);
xlabel('Period','FontSize',20);
ylabel('Market share','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_share0.pdf')


% C2 concentration index
figure('visible','off')
plot(c2_0,'-','LineWidth',2);
hold on
plot(c2_1,'r--','LineWidth',2);
xlabel('Period','FontSize',20);
ylabel('Market share','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_c2.pdf')

% Consumer surplus
figure('visible','off')
plot(csurp_0,'-','LineWidth',2);
hold on
plot(csurp_1,'r--','LineWidth',2);
xlabel('Period','FontSize',20);
ylabel('Surplus','FontSize',20);
legend('Baseline','Scenario 1','Location','northwest')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_csurp.pdf')

% Aggregate profits.
figure('visible','off')
plot(profits_0,'-','LineWidth',2);
hold on
plot(profits_1,'r--','LineWidth',2);
xlabel('Period','FontSize',20);
ylabel('Profits (US$ million)','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_profits.pdf')

% Average R&D expenditures.
figure('visible','off')
plot(avrd_0,'-','LineWidth',2);
hold on
plot(avrd_1,'r--','LineWidth',2);
xlabel('Period','FontSize',20);
ylabel('R&D expenditures (US$ million)','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_avRD.pdf')

% Average marginal costs.
figure('visible','off')
plot(avmgc_0,'-','LineWidth',2);
hold on
plot(avmgc_1,'r--','LineWidth',2);
xlabel('Period','FontSize',20);
ylabel('US$/watt','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_avmgc.pdf')

close all

table = [p0_0 share0_0 avp_0 avrd_0 c2_0 csurp_0 profits_0 avmgc_0;... 
         zeros(2,8);...
         p0_1 share0_1 avp_1 avrd_1 c2_1 csurp_1 profits_1 avmgc_1];

csvwrite('results/counter_output.csv',table);

end
%}
