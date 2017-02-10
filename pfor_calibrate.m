function [gamma,p0grid,Q] = pfor_calibrate(data)

% =========================================================================
% 
% p0_calibrate
%
% Calibrates the hazard of a foreign price change by minimizing the 
% distance between simulated price trajectories and the data. 
% Inputs are the foreign price data vector, initial hazard (g1) and initial
% probability (g2).
% Output includes (fixed) grid 'pfor', transition matrix 'Q', and implied ergodic 
% distribution over price states.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2015.
% =========================================================================

% Generate foreign price process objects.
disp('Start foreign price calibration.');
disp(' ');

% Parameters.
T = 25;
S = 500;
jumps = 1000;
rng(20027);
drws1 = rand(jumps,S);
drws2 = rand(jumps,S);

% Price grid.
p0grid = [
10.0450;        % avg 1999-2000.
7.1060;         % avg 2001-2005.
2.070;          % avg 2006-2009.
1.1244;         % avg 2010-2013.
0.75];        % proj 2014.

% Initial parameters: hazard, probability.
g10 = 0.4;
g11 = 0.4;
gstart = [g10;g11];
g = fminsearch(@(g) p0_dist(g,data,p0grid,T,S,jumps,drws1,drws2),gstart);
%g = fmincon(@(g) p0_dist(g,data,p0grid,T,S,jumps,drws1,drws2),[g10;g11],[],[],[],[],[0.00001;0.00001],[2;1]);

% Get output (moments+Q) given solved parmeters.
[d,mm,Q] = p0_dist(g,data,p0grid,T,S,jumps,drws1,drws2);

% Get moment from data.
[~,~,sd,mu,~,~] = ou_calib(data);

% Compute ergodic distribution.
[~,D,W] = eig(Q,'vector');
[~,i] = max(D);
erg = W(:,i)/sum(W(:,i));

disp(['Hazard rate: ', num2str(g(1))]);
disp(['Probability: ', num2str(g(2))]);
disp(['SD: ', num2str([sd mm(1)])]);
disp(['Mu: ', num2str([mu mm(2)])]);
disp(['Distance: ', num2str(d)]);

% Simulate price trajectories.
gamma = g(1);
p0t = zeros(T,S);

for s=1:S 
   p0t(:,s) = pfortraj(T,jumps,gamma,Q,p0grid,drws1(:,s),drws2(:,s)); 
end
p0t = mean(p0t,2);

save pfor_calibration p0grid Q erg g p0t

% Plot simulation vs data.
t = size(data,1);
data = [data;NaN(T-t,1)];
tt = (1999:1:1999+T-1);

figure%('visible','off')
plot(tt,data,'-','LineWidth',3.5);
hold on
plot(tt,p0t,'r--','LineWidth',3.5);
xlabel('Period','FontSize',20);
legend('Data','Average over S simulations','Location','northeast')
ylabel('$/Watt','FontSize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/p0_data-sim.pdf')

end


%% Compute distance metric.
function [d,model,Q] = p0_dist(g,data,pgrid,T,S,jumps,drws1,drws2)

% =========================================================================
% 
% p0_dist
%
% Simulates a set of trajectories for foreign price and computes the
% distance of the average over simulations and the data. 
%
% Written by Bernardo Diaz de Astarloa @ PSU 2015.
% =========================================================================


% Transition matrix.
% Build transition matrix.
Q = zeros(5,5);
Q(1,2) = 1;                 
Q(2,1) = 0.25;           % up from #4
Q(2,3) = 1 - Q(2,1);
Q(3,2) = 0.05;           % up from #3
Q(3,4) = 1 - Q(3,2);
Q(4,3) = g(2);           % up from #2
Q(4,5) = 1 - Q(4,3);
Q(5,4) = 1;              % up from #1

gamma = g(1);

%pt  = zeros(T,S);
sds = zeros(1,S);
mus = zeros(1,S);
for s=1:S
   p = pfortraj(T,jumps,gamma,Q,pgrid,drws1(:,s),drws2(:,s));
   [~,~,sd,mu,~,~] = ou_calib(p);
   sds(s)  = sd;
   mus(s)  = mu;
   %pt(:,s) = p; 
end

%pt = mean(pt,2);
%target = data;
%t = size(data,1);
%pt = pt(1:t);
%error = target - pt;

sd = mean(sds,2);
mu = mean(mus,2);
model = [sd;mu];
target = zeros(2,1);
[~,~,target(1),target(2),~,~] = ou_calib(data);
error = target - model;


d = norm(error)/norm(target);

end

%% Estimate AR(1) and match to OU process parameters.
function [a,b,sde,mu,sig,lambda] = ou_calib(data)

% Time-period: 1 if annual data.
delta = 1;

n = length(data)-1; 
Sx  = sum( data(1:end-1) );
Sy  = sum( data(2:end) );
Sxx = sum( data(1:end-1).^2 );
Sxy = sum( data(1:end-1).*data(2:end) );
Syy = sum( data(2:end).^2 );

% AR1 parameters: y(t+1) = ay(t) + b + eps.
a   = (n*Sxy - Sx*Sy)/(n*Sxx - Sx^2);
b   = (Sy - a*Sx)/n;
sde = sqrt((n*Syy - Sy^2 - a*(n*Sxy - Sx*Sy))/(n*(n-2)));

% Map to OU parameters.
lambda = -log(a)/delta;
mu  = b/(1-a);                                  % long run mean.
sig = sde*sqrt( (-2*log(a))/(delta*(1-a^2)) );

end

%% Compute foreign price trajectory.
function pt = pfortraj(T,jumps,gamma,Q,pfor,drw1,drw2) 

% =========================================================================
% 
% pfortraj
%
% Given uniform random draws for arrival times and parameters of the
% process, simulates a realization of the foreign price process and
% discretizes the series for T periods.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2014.
% =========================================================================

times_c = zeros(jumps,1);
p0_c    = zeros(jumps,1);

% Initial condition.
n     = 1;
p0    = 1;
p0_c(n)  = pfor(p0);

t = - log(drw1(n))/gamma;
n = n+1;

while t>0 && t<T 
    times_c(n) = t;
    ptrans = Q(p0,:);
    edge = [zeros(1,1) cumsum(ptrans,2)];
	draw = histc(drw2(n),edge);
    draw(end) = [];
    draw = find(draw==1);
    p0 =draw;
    
    % Update.
    p0_c(n) = pfor(p0);
    n = n+1;
    t = t - log(drw1(n))/gamma;
end

times_c = times_c(1:n-1,:);
p0_c = p0_c(1:n-1,:);


% Discretize continuous time data simulations.
pt = zeros(T,1);
t_lag = 0;

for t = 1:T
    t_ind = find(times_c(:,1)<t,1,'last');
    if isempty(t_ind)
        t_ind = t_lag;
    end
    if t_ind == t_lag
        if t==1
            pt(t) = p0_c(1);
        else
            pt(t) = pt(t-1);
        end
    else
        % There was a jump.
        % Integrate variables over subperiods.
        if t==1
            weight = [times_c(t_lag+2:t_ind);t] - times_c(t_lag+1:t_ind);
            pt(t)  = weight'*p0_c(t_lag+1:t_ind);
        else
            weight = [times_c(t_lag+1:t_ind);t] - [t-1;times_c(t_lag+1:t_ind)];
            pt(t)  = weight'*p0_c(t_lag:t_ind);
        end
    end
    
    % Update time period tracker.
    t_lag = t_ind;
end

end
