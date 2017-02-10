function avdata = averages(par,p1,x1,y1)

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

p0sim    = zeros(S,T);
avp      = zeros(S,T);
firms    = zeros(S,T);

drawany  = draws.any;
drawind  = draws.ind;
drawout  = draws.out;


for sim=1:S
    drw1 = drawany(:,sim);
    drw2 = drawind(:,sim);
    drw3 = drawout(:,sim);
    % Simulate industry: arrival times and jumps.
    outc = sim_cont(par,p1,x1,y1,di_0,drw1,drw2,drw3,T);

    % Fit to discrete time grid (years).
    % AD: aggregate data,
    % [ff,entry,exit,cs,ps,scrap,ecost,p0,q0,sh0,d,w,avprice,c2]
    % PD: panel data.
    [agg,~] = sim_disc(par,outc,T);
    p0sim(sim,:) = agg(:,8)';
    avp(sim,:)   = agg(:,13)';
    firms(sim,:) = agg(:,1)';
end

p0sim = nanmean(p0sim,1);
avp   = nanmean(avp,1);
firms = nanmean(firms,1);

avdata = [p0sim' avp' firms'];