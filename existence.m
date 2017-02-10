function [w_dvown,w_dvoth] = existence(par,V1,p1)

% =========================================================================
% 
% existence
%
% Check FOCs, second order condition of pricing FOC to establish existence,
% and non-decreasing/increasing properties of the value function.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2015.
% =========================================================================

% Call MEX function. 'checks' gives [foc,ex,dvown,dvoth].
[foc,ex,dvown,dvoth]=checks(par,V1,p1);
tol_foc = 1e-4;
norm=max(max(abs(foc)));

if norm>tol_foc;
    display(['WARNING: FOCs not satisfied to desired precision (',num2str(norm),').']);
else
    display('ALL GOOD: FOCs satisfied to desired precision');
end
if find(ex<0,1)>0;
    display('WARNING: Existence conditions not satisfied.');
else
    display('ALL GOOD: Existence conditions satisfied.');
end

w_dvoth = [];
w_dvown = [];

if find(dvown<0,1)>0;
    display('WARNING: Value function is decreasing in own w at some states.');
    [~,J] = find(dvown<0);
    J = unique(J);
    w_dvown = par.state(:,J);
end
if find(dvoth>0,1)>0;
    display('NOTE: Value function is increasing in others w at some states.');
    [~,J] = find(dvoth>0);
    J = unique(J);
    w_dvoth = par.state(:,J);
end


% MATLAB version.
%{
function [foc,check_foc,check_ex,check_dv,states_dv] = existence(p1,V1,par)

foc = zeros(par.N,par.S);   
check_foc = zeros(par.N,par.S);             
check_ex = zeros(par.N,par.S);             
check_dvown = zeros(par.N,par.S);
check_dvoth = zeros(par.N,par.S);
states_dv = [];
a = par.alpha;
p0size = par.p0size;
p0vec  = par.pfor;
mkt    = par.mkt;
N = par.N;
D = par.D;
M = par.M;
state = double(par.state);
binom = double(par.binom);

for w = 1:par.S
    % Load state.
    di = state(:,w);
    nstar = sum(di(2:N+1)~=M+1);
    p   = p1(:,w);
    v   = V1(:,w);
    d0  = di(1);
    d   = ceil(d0/p0size);
    p0  = d0 - (d-1)*p0size;
    d   = mkt(d);
    p0  = p0vec(p0);
    c   = par.mgc(di(2:nstar+1));
    
    % Shares, quantities and hazard derivatives.
    den = exp(-a*(p(1:nstar) - p0));
    s   = exp(-a*(p(1:nstar) - p0))./(1 + sum(den)); % shares
    q   = d*s;
    hazp  = par.eta1*par.eta2*q.^(par.eta2-1);
    hazpp = par.eta1*par.eta2*(par.eta2-1)*q.^(par.eta2-2);

    for j = 1:nstar
        % Compute DeltaVji for all i.
        vj = v(j);
        Dv = zeros(nstar,1);
        for i = 1:nstar
            dip = di;
            dip(i+1) = min(dip(i+1)+1,M);
            [ww,pos] = encode(dip,binom,N,D);
            vup = V1(pos,ww);
            Dv(i) = vup(j) - vj;
        end
        pj = p(j);
        qj = q(j);
        cj = c(j);
        sj = s(j);
        si = s;
        si(j) = 0;
        diffj = Dv(j);
        diffi = Dv;
        diffi(j) = 0;
        hazpj  = hazp(j);
        hazpi = hazp;
        hazpi(j) = 0;
        hazppj = hazpp(j);
        hazppi = hazpp;
        hazppi(j) = 0;
        
        % Check first order condition.
        f = 1/a - (1-sj)*(pj-cj) - (1-sj)*hazpj*diffj + (hazpi.*si)'*diffi;
        foc(j,w) = f;
        check_foc(j,w) = abs(f)<1e-4;
        
        % Check second order condition.
        lhs = 1/(a*qj);
        rhs = hazppj*(1-sj)^2*diffj + (hazppi.*(si.^2))'*diffi;
        check_ex(j,w) = lhs>rhs;
        
        % Check monotonicity of value function.
        check_dvown(j,w) = diffj>=0;             % own effect.
        check_dvoth(j,w) = sum(diffi>0)==0;      % others effects.
    end
    if nstar<N
        % Fill in non-active firms with 2s not to confound with checks.
        check_foc(nstar+1:end,w) = 2*ones(N-nstar,1);
        check_ex(nstar+1:end,w) = 2*ones(N-nstar,1);
        check_dvown(nstar+1:end,w) = 2*ones(N-nstar,1);
        check_dvoth(nstar+1:end,w) = 2*ones(N-nstar,1);
    end
end
%}
