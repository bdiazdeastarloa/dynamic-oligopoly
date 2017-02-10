s2 = zeros(par.S,1);
d0 = 15;
i=1;
for s = 1:par.S
    if par.state(4,s)==M+1 && par.state(3,s)<M+1 && par.state(1,s)==d0
        s2(i) = s;
        i = i+1;
    end
end
s2 = s2((s2>0));
ss = par.state(2:3,s2);
pp = p1(1:2,s2);

psurf  = zeros(M,M);
jj = [1;0];
for i=1:M                     % firm 1's state
    for j=1:M                 % firm 2's state
        nw=[i;j];
        temp = sortrows([nw,jj]);   % sort decreasing  
        d = temp(:,1);
        [~,l] = max(temp(:,2));                       
        for s=1:size(ss,2)
            if ss(1,s)==d(1) && ss(2,s)==d(2);
                psurf(i,j) = pp(l,s);
            end
        end
    end
end

figure(2)
mesh(psurf);