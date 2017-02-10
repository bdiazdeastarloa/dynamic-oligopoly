%% Firm trajectories

T = 11;                                      % # of periods to simulate
S = 500;                                     % # of simulations to run
M = par.M;
p0size = par.p0size;

mkt0 = 1;                                    % initial demand state
p0   = 1;                               % initial foreign price state
w0   = [1 2 3 3 3 M+1 M+1 M+1]';             % initial industry structure
d0   = p0+(mkt0-1)*p0size;                   % index of exogenous shock
di_0 = [d0;w0];

savedraws = sprintf('draws_T%dS%d.mat',[T,S]); 
load(savedraws);

% Pick vector or random draws;
sim = 182;                 
drw1 = draws.any(:,sim);
drw2 = draws.ind(:,sim);
drw3 = draws.out(:,sim);

% Continuous time simulation.
outc = sim_cont(par,p1,x1,y1,di_0,drw1,drw2,drw3,T);
% Discretize simulation: only firm-level outcomes.
[agg,panel] = sim_disc(par,outc,T);

% Preallocate (there are nf firms).
fprices = NaN(T,outc.nf);
fmgc    = NaN(T,outc.nf);
fquant  = NaN(T,outc.nf);
fprof   = NaN(T,outc.nf);
for i=1:outc.nf
    in   = panel(panel(:,1)==i,2);
    temp = panel(panel(:,1)==i,4);                               
    fprices(in(:,1),i) = temp;
    temp = panel(panel(:,1)==i,9);
    fmgc(in(:,1),i) = temp;
    temp = panel(panel(:,1)==i,6);
    fquant(in(:,1),i) = temp;
    temp = panel(panel(:,1)==i,8);      % profits
    temp2 = panel(panel(:,1)==i,11);    % entry costs
    temp3 = panel(panel(:,1)==i,13);    % rd expenditures
    fprof(in(:,1),i) = temp - temp3;        
end


%% Average selling price.
figure('visible','off')
plot(fprices(:,1),'-','LineWidth',3.5);
hold on
plot(fprices(:,2),'--','LineWidth',2);
hold on
plot(fprices(:,3),'r-','LineWidth',2);
hold on
plot(fprices(:,4),'r--','LineWidth',2);
hold on
plot(fprices(:,5),'b-','LineWidth',2);
hold on
plot(fprices(:,6),'b:','LineWidth',3);
hold on
plot(fprices(:,7),'g-.','LineWidth',2);
%hold on
%plot(fprices(:,8),'g-','LineWidth',3);
%hold on
%plot(fprices(:,9),'g--','LineWidth',3);
xlabel('Period','FontSize',20);
ylabel('$/watt','FontSize',20);
%legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7','Firm 8','Firm 9','Location','northeast')
legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/ftraj_prices.pdf')


%% Marginal cost.
figure('visible','off')
plot(fmgc(:,1),'-','LineWidth',3.5);
hold on
plot(fmgc(:,2),'--','LineWidth',2);
hold on
plot(fmgc(:,3),'r-','LineWidth',2);
hold on
plot(fmgc(:,4),'r--','LineWidth',2);
hold on
plot(fmgc(:,5),'b-','LineWidth',2);
hold on
plot(fmgc(:,6),'b:','LineWidth',3);
hold on
plot(fmgc(:,7),'g-.','LineWidth',2);
%hold on
%plot(fmgc(:,8),'g-','LineWidth',3);
%hold on
%plot(fmgc(:,9),'g--','LineWidth',3);
xlabel('Period','FontSize',20);
ylabel('$/watt','FontSize',20);
%legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7','Firm 8','Firm 9','Location','northeast')
legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/ftraj_mgc.pdf')


%% Quantities.
figure('visible','off')
plot(fquant(:,1),'-','LineWidth',3.5);
hold on
plot(fquant(:,2),'--','LineWidth',2);
hold on
plot(fquant(:,3),'r-','LineWidth',2);
hold on
plot(fquant(:,4),'r--','LineWidth',2);
hold on
plot(fquant(:,5),'b-','LineWidth',2);
hold on
plot(fquant(:,6),'b:','LineWidth',3);
hold on
plot(fquant(:,7),'g-.','LineWidth',2);
%hold on
%plot(fquant(:,8),'g-','LineWidth',3);
%hold on
%plot(fquant(:,9),'g--','LineWidth',3);
xlabel('Period','FontSize',20);
ylabel('MW','FontSize',20);
%legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7','Firm 8','Firm 9','Location','northeast')
legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/ftraj_ship.pdf')


%% Profits.
figure('visible','off')
plot(fprof(:,1),'-','LineWidth',3.5);
hold on
plot(fprof(:,2),'--','LineWidth',2);
hold on
plot(fprof(:,3),'r-','LineWidth',2);
hold on
plot(fprof(:,4),'r--','LineWidth',2);
hold on
plot(fprof(:,5),'b-','LineWidth',2);
hold on
plot(fprof(:,6),'b:','LineWidth',3);
hold on
plot(fprof(:,7),'g-.','LineWidth',2);
%hold on
%plot(fprof(:,8),'g-','LineWidth',3);
%hold on
%plot(fprof(:,9),'g--','LineWidth',3);
xlabel('Period','FontSize',20);
ylabel('US$ million','FontSize',20);
%legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7','Firm 8','Firm 9','Location','northeast')
legend('Firm 1','Firm 2','Firm 3','Firm 4','Firm 5','Firm 6','Firm 7')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/ftraj_profits.pdf')


%% Foreign price.
figure('visible','off');
[xx,price0,share0] = plotyy(1:T,agg(:,9),1:T,agg(:,11),'plot');
set(get(xx(1),'Ylabel'),'String','Foreign price (US$/watt)','Fontsize',20) 
set(get(xx(2),'Ylabel'),'String','Imports market share','Fontsize',20)
set(xx(2),'xticklab',[],'xtick',[])
set(xx(2),'Fontsize',20)
set(xx(1),'Fontsize',20)
set(price0,'LineStyle','-','Color','b','Linewidth',2.5)
set(share0,'LineStyle','--','Color','r','Linewidth',2.5)
xlabel('Period','Fontsize',20)
legend('Foreign price','Market share','Location','northeast');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/traj_forprice.pdf')
