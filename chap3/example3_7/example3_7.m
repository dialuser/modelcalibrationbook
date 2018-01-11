%{
Copyright (C) {{2015}}  {{
%MODEL CALIBRATION AND PARAMETER ESTIMATION,  N.Z. Sun and A.Y. Sun}}

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%!!! Tested for Matlab R2011
%}

%Author: Alex
%Date: $12242013$
%Purpose: HyMod driver
%All units are in [mm]
%For Example 3.7
%Dependency:
%This example requires Matlab Global Optimization Toolbox 
%==========================================================================
clear all;
global dailyPrecip dailyPotEvapTrans dailySS;
global t0 tEnd;
%% User-defined Parameter
%Maximum combined contents of all stores, [L]
par.cmax  = [150 ,350];
%Scaled distribution function shape parameter [-] 
par.bexp  = [0.10, 1.5];
%Quick/slow split parameter   [-]
par.fQuickFlow = [0.60,0.99];
%Slowflow routing tank's rate parameter or residence time [days]
par.Rs    = [0.01, 0.1];
%Quickflow routing tanks' rate parameter or residence time [days]
par.Rq    = [0.20,0.7];

%load input data
%This is area of Guadalupe basin, not used
basinArea = 3406*1e6; %[m^2]

fid = fopen('leaf_rv.in', 'r');
%jump 35 line
for i=1:35
    fgets(fid);
end
data=textscan(fid, '%8c%10n%10n%10n%10n%10n');
% yy = data{1}(:);
% mm = data{2}(:);
% dd = data{3}(:);
% ssdates = datenum(yy, mm, dd); 

dailyPrecip = data{2}(:); %[mm/day]
dailyPotEvapTrans = data{3}(:); %[mm/day]
dailySS = data{4}(:); %in [mm/day]
%convert to m3/day if necessary
iconvert = 0;
if (iconvert)
    dailyPotEvapTrans = dailyPotEvapTrans./1000*basinArea;
    dailyPrecip = dailyPrecip./1000*basinArea;
    dailySS = dailySS./1000*basinArea; % [m^3/day]
end
fclose(fid);

%calibrate SS using data from 1948/1/1 to 12/31/1979
nRec = length(dailySS);
trainRange = 1:(11723-35+1);
valRange = (trainRange(end)+1):(19394-35);
t0 = trainRange(1);
tEnd = trainRange(end);
[testval, outflow] = hymod([300, 1.0, 0.5, 0.5, 0.5]);
figure(1),clf;
hold on;
plot(outflow, 'r-', 'LineWidth', 1.5);
plot(dailySS(trainRange), 'b--');
hold off;

%% set up GA
options = gaoptimset(...
                     'MutationFcn',@mutationadaptfeasible, ...
                     'Display', 'diagnose', ...
                     'TolFun', 1e-7, ...
                     'CrossoverFraction', 0.8, ...
                     'SelectionFcn', @selectiontournament, ...
                     'PlotFcns', @gaplotbestf);

rng = RandStream.create('mt19937ar','seed',51899);
RandStream.setDefaultStream(rng);

% Define parameer bounds
lb = [par.cmax(1) par.bexp(1) par.fQuickFlow(1) par.Rs(1) par.Rq(1)];
ub = [par.cmax(2) par.bexp(2) par.fQuickFlow(2) par.Rs(2) par.Rq(2)];

[thetaparam,fval,exitflag,output] = ga(@hymod,5,[],[],[],[],...
                                    lb,ub,[], options);
save('ga_multi.mat', 'thetaparam');                                
% Do training
nSol = size(thetaparam,1);
trainObj = zeros(nSol,2);
valObj = zeros(nSol,2);

for i=1:nSol
    [trainObj(i,:), outflow] = hymod(thetaparam(i,:));
end
figure(1), clf;
subplot(1,2,1);
hold on;
plot(dailySS(trainRange), 'o');
plot(outflow, 'b-', 'LineWidth', 1.5);
hold off;

% Do validation
t0 = valRange(1);
tEnd = valRange(end);
for i=1:nSol
    [valObj(i,:), outflowVal] = hymod(thetaparam(i,:));
end
meanobs = mean(dailySS(valRange));
nsc = 1 - sum((outflowVal-dailySS(valRange)).^2)/sum((dailySS(valRange)-meanobs).^2);
subplot(1,2,2);
hold on;
plot(dailySS(valRange), 'o');
plot(outflowVal, 'b-', 'LineWidth', 1.5);
hold off;
hold off;