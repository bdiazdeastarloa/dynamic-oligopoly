function [np,nv] = init(par)

state = double(par.state);
agg   = par.agg;
alpha = par.alpha;
rho   = par.rho;
mkt   = par.mkt;
mgc   = par.mgc;
pfor  = par.pfor;
N     = par.N;
S     = par.S;
M     = par.M;

usedogleg = 1;
np = zeros(N,S);
nv = zeros(N,S);

for s = 1:S
    
    % Obtain unique state.
    di = state(:,s);
    
    % Obtain number of active firms.
    Nstar = sum(di(2:N+1)~=M+1);
    
    if Nstar>0
        c  = mgc(di(2:Nstar+1));
        p0 = pfor(agg(di(1),2));
        % Check we are working with (Nx1) vectors.
        if size(c,2)~=1
            c = c';
        end
        pars.c  = c;
        pars.p0 = p0;
        pars.a  = alpha;
        if usedogleg
            MaxIter = 400; TolFun = 1e-6; TolX = 1e-6; 
            options = [1,TolFun,TolX,TolFun,MaxIter,1e-6];
            [pstar,info] = SDogLeg('foc_static',pars,c,options);
            if info(6)>3
                warning('There is a problem with SDogLeg!');
                pause
            end
        else
            options = optimset('Display','iter');
            [pstar,~,exitflag] = fsolve('foc_static',c,options,pars);
            if exitflag<=0
                warning('There is a problem with fsolve!');
                pause
            end
        end
        
        % Compute profits.
        d = exp(-alpha*(pstar-p0));               
        sh = d./(1+sum(d));
        m = mkt(agg(di(1),1));
        pr = m*sh.*(pstar - c);
        
        % Update value and price.
        np(1:Nstar,s) = pstar;
        nv(1:Nstar,s) = pr./rho;
    end
    
    % if (mod(s,5000) == 0)
    %     disp(['Solving state ' int2str(s) '.']);
    % end        
end