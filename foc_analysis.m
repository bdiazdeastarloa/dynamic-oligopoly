function foc_analysis
close all

p1  = -100:1:500;
p2  = [5.3457;5.3457];
pp = struct;
pp.p0 = 3.8891;
pp.alpha = 2.779;
pp.eta = 0.0009;
pp.gam = 0.004;
pp.nu  = 0.6;
pp.cj = 3.8891;
pp.M = 160;
deltav = 1e03*[1.127919864830376; -0.002800336000167;  -0.002800336000167];

i = 1;
y1  = zeros(i,size(p1,2));
y2  = zeros(i,size(p1,2));
y3  = zeros(i,size(p1,2));

for j=1:i;
    for in=1:size(p1,2)
        y1(j,in) = foc(p1(in),p2,pp,deltav(:,j),1);
        y2(j,in) = foc(p1(in),p2,pp,deltav(:,j),2);
        y3(j,in) = foc(p1(in),p2,pp,deltav(:,j),3);
    end
end

%xx = startval(3.8900,p2,pp,deltav(:,j));

figure(1)
plot(p1,y1(1,:));
hold on
plot(p1,y2(1,:),'r');
plot(p1,y3(1,:),'g');
axis([-10 100 -5 100]);

end

function f = foc(p1,p2,pp,deltav,h)

    p     = [p1;p2];
    alpha = pp.alpha;
    eta   = pp.eta;
    gam   = pp.gam;
    nu    = pp.nu;
    p0    = pp.p0;
    
    denom = 1+sum(exp(- alpha*(p-p0)));               
    s     = exp(- alpha*(p-p0))./denom;
    cj   = pp.cj;
    q    = pp.M*s;
    
    if h==1
        hq = eta;             % derivative of hazard function
    elseif h==2
        hq = gam./((1+gam*q).^2);
    else
        hq = eta*nu*(q.^(nu-1));
    end
    hqs = hq.*s;
    
    f     = 1 - alpha*((1-s(1))*(p1-cj) - hq(1)*deltav(1) + hqs'*deltav);
    
end

function z = startval(p1,p2,pp,deltav,h)

    p     = [p1;p2];
    alpha = pp.alpha;
    eta   = pp.eta;
    nu    = 1;
    p0    = pp.p0;
    
    denom = 1+sum(exp(- alpha*(p-p0)));               
    s     = exp(- alpha*(p-p0))./denom;
    cj   = pp.cj;
    q    = pp.M*s;

    q1    = eta*nu*q.^(nu-1);             % derivative of hazard function
    q2    = q1.*s;

    z     = cj+(q2'*deltav - 1/alpha)/(1-s(1));
    
end

