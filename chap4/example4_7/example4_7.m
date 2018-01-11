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

function example4_7()
%Example 4_7
%Author: Alex Sun
%Date: $20140102$
%Purpose: Estimate HyMOD parameters using MCMC
%Change isimu = 1 to regenerate all results
%Dependency:
% Requires mcmcstat (included)
%=====================================================================
clear all;
global dailyPrecip dailyPotEvapTrans dailySS;
global t0 tEnd;
%% User-defined Parameter
%Maximum combined contents of all stores, [mm]
par.cmax  = [150, 350];
%Scaled distribution function shape parameter [-] 
par.bexp  = [0.10, 1.5];
%Quick/slow split parameter   [-]
par.fQuickFlow = [0.6,0.99];
%Slowflow routing tank's rate parameter or residence time [days]
par.Rs    = [0.01, 0.1];
%Quickflow routing tanks' rate parameter or residence time [days]
par.Rq    = [0.20,0.7];

fid = fopen('leaf_rv.in', 'r');
%number of lines to jump
headers=17;
for i=1:headers
    fgets(fid);
end
data=textscan(fid, '%8c%10n%10n%10n%10n%10n');
dailyPrecip = data{2}(:); %[mm/day]
dailyPotEvapTrans = data{3}(:); %[mm/day]
dailySS = data{4}(:); %in [mm/day]
fclose(fid);

%calibrate SS using data from 1952/1/1 to 12/31/1962
nRec = length(dailySS);
trainRange = 1478-headers+1:(11723-headers+1);
valRange = (trainRange(end)+1):(19394-headers+1);
%set start and end time for calibration
t0 = trainRange(1);
tEnd = trainRange(end);
% Define parameer bounds
lb = [par.cmax(1) par.bexp(1) par.fQuickFlow(1) par.Rs(1) par.Rq(1)];
ub = [par.cmax(2) par.bexp(2) par.fQuickFlow(2) par.Rs(2) par.Rq(2)];

%% set up MCMC
%x is used to deonte parameter to be estimated, 
%x = (cmax, bexp, fQuickflow, Rs, Rq)
%Initial guess
nDim = 5; %dimension of the unknown
x0 =  [220 0.5 0.9 0.02 0.6]; %initial guess
mu0 = [220 0.5 0.9 0.02 0.6]; %mean of prior
sigma0 =  [100 1.0 0.1 0.1 0.3]; %std dev of prior
sigmaq =  [100 1.0 0.1 0.1 0.3]; %std dev of proposal
covmat = zeros(nDim);
for i=1:nDim
    covmat(i,i) = sigmaq(i).^2;
end

addpath([pwd '/mcmcstat']);

%param structure, {'name', initial, min, max, prior_mean, prior_sigma}
params = {
    {'cpar', x0(1), lb(1), ub(1),  mu0(1), sigma0(1)}
    {'bexp', x0(2), lb(2), ub(2),  mu0(2), sigma0(2)}
    {'alpha', x0(3), lb(3), ub(3), mu0(3), sigma0(3)}
    {'Rs', x0(4), lb(4), ub(4),    mu0(4), sigma0(4)}
    {'Rq', x0(5), lb(5), ub(5),    mu0(5), sigma0(5)}
};
data={};
data.xdata = dailySS(trainRange);

model.ssfun  = @mylikelihood;
model.sigma2 = 0.01^2;
model.N =  length(trainRange);

options.nsimu = 5500;
options.updatesigma = 1;
options.qcov = covmat;
options.burnintime = 500;
options.waitbar = 1;
%set to 1 to regenerate, set to 0 to load existing results
isimu = 1;
if (isimu)
    rng = RandStream.create('mt19937ar','seed',23999);
    RandStream.setGlobalStream(rng);

    [res,chain,s2chain] = mcmcrun(model,data,params,options);
    save('example4_7.mat', 'res', 'chain', 's2chain', 'model', ...
        'data', 'params', 'options');
else
    load('example4_7.mat');
end
postmean = res.mean;

figure(1), clf;

subplot(2,3,1);
chain = chain(options.burnintime:end,:);
rhist(chain(:,1),20); title('(a)');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8],'EdgeColor','k');
xlabel('C_{max}');
ylabel('Frequency');

subplot(2,3,2);
rhist(chain(:,2),20), title('(b)');
xlabel('b');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8],'EdgeColor','k');
ylabel('Frequency');

subplot(2,3,3);
rhist(chain(:,3),20), title('(c)');
xlabel('\alpha');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8],'EdgeColor','k');
ylabel('Frequency');

subplot(2,3,4);
rhist(chain(:,5),20), title('(d)');
xlabel('R_q');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8],'EdgeColor','k');
ylabel('Frequency');

subplot(2,3,5);
rhist(chain(:,4),20), title('(e)');
xlabel('R_s');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8],'EdgeColor','k');
ylabel('Frequency');


figure(2), clf;
disp(postmean);
subplot(1,5,1);
plot(chain(:,1));
subplot(1,5,2);
plot(chain(:,2));
subplot(1,5,3);
plot(chain(:,3));
subplot(1,5,4);
plot(chain(:,5));
subplot(1,5,5);
plot(chain(:,4));


%% This is the likelihood
function [res] =mylikelihood(x, data)
        %call hymod to calculate objective function
        res = hymod(x);
end
end