function count_out

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
% 
% Output of counterfactual analysis. Builds figures and tables with time
% series.
%
% Needs counterfactual simulations to have been generated first (in
% 'counterfactuals.m'.
% Output:
% - Figures saved as PDFs. 
% - Table with output variables saved as an XLS file.
%
% Written by Bernardo Diaz de Astarloa @ PSU March 2015
% =========================================================================

load results/simcount_c0.mat
nfirms_0  = simcount_m.nfirms;
exit_0    = simcount_m.exit; 
entry_0   = simcount_m.entry;
avp_0     = simcount_m.avprice;
p0_0      = simcount_m.pfor;
share0_0  = simcount_m.share0;
omega_0   = simcount_m.omega;
c2_0      = simcount_m.c2;
csurp_0   = simcount_m.csurplus;
profits_0 = simcount_m.netprofits;
avrd_0    = simcount_m.avrd;
rdjump_0  = simcount_m.rdjump;
lbdjump_0 = simcount_m.lbdjump;
d_cs_0    = simcount_m.d_cs;
d_psx_0   = simcount_m.d_psx;
d_sc_0    = simcount_m.d_sc;
d_ec_0    = simcount_m.d_ec;

load results/simcount_c1.mat
nfirms_1  = simcount_m.nfirms;
exit_1    = simcount_m.exit; 
entry_1   = simcount_m.entry;
avp_1     = simcount_m.avprice;
p0_1      = simcount_m.pfor;
share0_1  = simcount_m.share0;
omega_1   = simcount_m.omega;
c2_1      = simcount_m.c2;
csurp_1   = simcount_m.csurplus;
profits_1 = simcount_m.netprofits;
avrd_1    = simcount_m.avrd;
rdjump_1  = simcount_m.rdjump;
lbdjump_1 = simcount_m.lbdjump;
d_cs_1    = simcount_m.d_cs;
d_psx_1   = simcount_m.d_psx;
d_sc_1    = simcount_m.d_sc;
d_ec_1    = simcount_m.d_ec;

% Domestic average selling price
figure('visible','off')
plot(avp_0,'-','LineWidth',3.5);
hold on
plot(avp_1,'r--','LineWidth',3.5);
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
plot(p0_0,'-','LineWidth',3.5);
hold on
plot(p0_1,'r--','LineWidth',3.5);
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
plot(share0_0,'-','LineWidth',3.5);
hold on
plot(share0_1,'r--','LineWidth',3.5);
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
plot(c2_0,'-','LineWidth',3.5);
hold on
plot(c2_1,'r--','LineWidth',3.5);
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
plot(csurp_0,'-','LineWidth',3.5);
hold on
plot(csurp_1,'r--','LineWidth',3.5);
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
plot(profits_0,'-','LineWidth',3.5);
hold on
plot(profits_1,'r--','LineWidth',3.5);
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
plot(avrd_0,'-','LineWidth',3.5);
hold on
plot(avrd_1,'r--','LineWidth',3.5);
xlabel('Period','FontSize',20);
ylabel('R&D expenditures (US$ million)','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_avRD.pdf')

% Productivity.
figure('visible','off')
plot(omega_0,'-','LineWidth',3.5);
hold on
plot(omega_1,'r--','LineWidth',3.5);
xlabel('Period','FontSize',20);
ylabel('Productivity','FontSize',20);
legend('Baseline','Scenario 1','Location','northeast')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/count1_omega.pdf')

close all

table = [p0_0 share0_0 avp_0 avrd_0 c2_0 csurp_0 profits_0 rdjump_0 lbdjump_0 nfirms_0 exit_0 entry_0;... 
         zeros(2,12);...
         p0_1 share0_1 avp_1 avrd_1 c2_1 csurp_1 profits_1 rdjump_1 lbdjump_1 nfirms_1 exit_1 entry_1];
welf  = [d_cs_0 d_psx_0 d_sc_0 d_ec_0;...
           d_cs_1 d_psx_1 d_sc_1 d_ec_1];
       
csvwrite('results/counter_output.csv',table);
csvwrite('results/counter_welfare.csv',welf);

end