function [D,W,error] = distance(param)

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
%
% distance:
% - assign current parameter guess.
% - call main routine.
% - assign data moments.
% - compute distance metric.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2014.
% =========================================================================

format long;

% Call main routine (solution + simulation).
% main uses Y, generates moments stored in MM.
parguess = param; %#ok
dtime = tic;

% Call main routine (solve for an equilibrium and compute moments).
main;

% Data.
p_pc90Data    = 3.464716;
p_ratioData   = 1.527128;            % Price 90th/5th percentile ratio.
enrateData    = 0.0757576;           % Mean (annual) entry rate.
exrateData    = 0.0454545;           % Mean (annual) exit rate.
p_acData      = 0.3811814;           % (Log) Price one-period autocorrelation.
revar1Data    = 0.7239576;           % Revenues AR1 regression (includes trend).
rd1Data       = 14.32406;            % R&D/Revenue for 1st firm in share increase.
rd2Data       = 15.09008;            % R&D/Revenue for 2nd firm in share increase.
rdar1Data     = 0.7711728;           % R&D AR1 regression (includes trend).
prod_qrd1Data = 0.0098201;           % Log price on acc. shipments + acc. RD, shipments.
prod_qrd2Data = -0.067045;           % Log price on acc. shipments + acc. RD, RD.

% List of targets to match. 
% (NOTE: order of vars matters as it is the same as bootcov below.)
Mlist = [
p_pc90Data;
p_ratioData;
enrateData;
exrateData;
p_acData;
rdar1Data;
revar1Data;
prod_qrd1Data;
prod_qrd2Data;           
rd1Data;
rd2Data];

% Bootstrapped standard errors and weight matrix.
load bootcov;
Wopt = inv(bootcov);            % Optimal weighting matrix.
W = diag(diag(Wopt));           % Actual weighting matrix.

% Compute loss function (MM is an output of 'main' above).
data = Mlist; 

error = data-MM;
D     = error'*W*error;
Duw   = norm(error)/norm(data);

nanflag = isnan(D); 
if nanflag>0;
    D = 1000; 
end

clear bootcov

labels = [
'p90    ';
'p90-5  ';
'En     ';
'Ex     ';
'p_AR1  ';
'RD_AR1 ';
'RV_AR1 ';
'LR_Q   ';
'LR_RD  ';
'RD_f1  ';
'RD_f2  '];

% Print Diagnostics
compare = cat(2,data,MM);

fprintf('==============================================================\n');
fprintf('\nParameters values:\n');
disp(   num2str(param));
fprintf('\nData - Model comparison:\n');
disp(   [labels num2str(compare)]);   
fprintf('\nLoss (weighted/unweighted):\n');
disp(   [D Duw]);
fprintf('This round time: %f\n', toc(dtime)/60);
fprintf('==============================================================\n');

