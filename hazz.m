%% Compute hazards of jumps.
function [hj,sh] = hazz(par,di,q,x,y)

% Computes 'hj', a vector of event-specific hazards, and 'sh', the hazard
% of a jump in the state of the industry due to any shock. 
% Note: hazard vectors are all column vectors.

in1 = (di(2:end)==par.M);               % can't increase productivity.
in2 = (di(2:end)==1);                   % can't decrease productivity.
in3 = find(di(2:end)<par.M+1,1,'last'); % # of active firms.
en  = find(di(2:end)==par.M+1,1);       % entrant location.

% LBD hazard.
h1 = par.eta1*q.^par.eta2;
% h1 = (par.eta1*q)./(10+par.eta1*q);          
h1(in1) = 0;
% R&D hazard.
h2 = par.alpha1*x(1:in3).^par.alpha2;
% h2 = par.alpha2*x(1:in3)./(1+par.alpha2*x(1:in3));  
h2(in1)= 0;
% Negative prod. shock hazard.
h3 = ones(in3,1)*par.delta;            
h3(in2) = 0;
% Exit hazard.
h4 = par.etax1*y(in3);   
% Entry hazard.
h5 = par.etae1*y(en);                  
% Aggregate shock hazard.
h6 = par.gamma;                              

if isempty(h1) 
    h1=zeros(in3,1);
end
if isempty(h2)
    h2=zeros(in3,1);
end
if isempty(h3)
    h3=zeros(in3,1);
end
if isempty(h4)
    h4=zeros(in3,1);
end
if isempty(h5)
    h5=0;
end

% h1 and h2 are not combined so that we can keep track of LBD vs R&D.
% Note the dim of hj is (4j+2,1).
hj = [h1;h2;h3;h4;h5;h6];
sh = sum(h1)+sum(h2)+sum(h3)+sum(h4)+h5+h6;

end
