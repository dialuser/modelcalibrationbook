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

%Dependencies: Requires Matlab ANN Toolbox

%Author: Alex Sun
%Example 8.4
%Purpose: Demonstrate ANN for daily streamflow forecasting
%Data: 
%Streamflow: USGS gauge 08033500, Neches Rv nr Rockland, TX
%http://waterdata.usgs.gov/nwis/nwisman/?site_no=08033500
%3,636 square miles
%Duration 1948/1/1-2003/12/31
%Other data: MOPEX 
%ftp://hydrology.nws.noaa.gov/pub/gcip/mopex/US_Data/Us_438_Daily/08055500.dly
%==========================================================================
%% Load daily streamflow 
%
clear all;
close all;
fid = fopen('08033500_Q.dat', 'r');
res = textscan(fid, '%s %f');
fclose(fid);
alldates=res{1}(:); %date strings
ss0 = res{2}(:);
% Transform the streamflow data using box-cox transform
[ss1, lambda] = boxcox(ss0);

%% Load daily precipitation, pet, tmax, and tmin data for the same gauge 
% (we don't read streamflow from this dataset because it is incomplete)
fid = fopen('08033500_daily.dat', 'r');
res = textscan(fid, '%f %f %*f %f %f');
fclose(fid);
precip = res{1}(:);
pet = res{2}(:);
tmax = res{3}(:);
tmin = res{4}(:);

%% Scaling data
nrec = length(ss0);
% linear scaling
XMIN = -0.99;
XMAX = 0.99;

[ss, SS] = mapminmax(ss1', XMIN, XMAX);    
tmax = mapminmax(tmax', XMIN, XMAX);
tmin = mapminmax(tmin', XMIN, XMAX);
precip = mapminmax(precip', XMIN, XMAX);
%% Uncomment to do correlation analysis
corranalysis = false;
if (corranalysis)
    figure(1), clf;
    subplot(3,2,1); 
    %autocorrelation
    autocorr(ss,60);title('Q autocorrelation');
    %cross-correlation
    subplot(3,2,2); 
    crosscorr(ss, precip); title('Q vs. precip');
    %
    subplot(3,2,3); 
    crosscorr(ss, tmax); title('Q vs. tmax');
    subplot(3,2,4); 
    crosscorr(ss, tmin); title('Q vs. tmin');
    subplot(3,2,5); 
    crosscorr(ss, tavg); title('Q vs. tavg');
end

%% Construct ANN for 1-day ahead prediction
%reset random seed
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',122391));
%
offset=10;
nday = 1;
allRange=offset:nrec-nday;
inputvec = [                
                precip(allRange);
                precip(allRange-1);
                precip(allRange-2);
                tmax(allRange); 
                tmax(allRange-1); 
                tmin(allRange);
                tmin(allRange-1); 
            ];
for iset= 0:offset-1
     inputvec=[inputvec;
                ss(allRange-iset)];
end
targetvec = ss(allRange+nday);

% Set network parameters

% Turn on testHiddenNeuron to test the effect of number of neurons
% on test data performance, from which the best num of hidden neurons
% is found
testHiddenNeuron = false;
if (testHiddenNeuron)
    neuronTests=2:10;
    neuronError=[];
    perf=[];
    for i = neuronTests
        net = newff(inputvec,targetvec, i, {'tansig'});
        %Do ANN ensemble if necessary, which is not done in this example
        %To do ANN ensemble, use net=init(net) to re-initiate net each time

        net.divideFcn='divideblock';
        net.divideParam.trainRatio = 0.65;  % Adjust as desired
        net.divideParam.valRatio = 0.15;  % Adjust as desired
        net.divideParam.testRatio = 0.20;  % Adjust as desired
        net.trainFcn = 'trainlm';
        % Train and Apply Network
        [net,tr] = train(net,inputvec,targetvec);
        ind = tr.testInd;
        P = inputvec;
        %simulate for testing period
        [~,~,~,~,perfval] = sim(net,P(:,ind));    
        neuronError=[neuronError perfval]    
    end
        
    AICtest = length(ind)*log(neuronError)+2*neuronTests ...
        +2.*neuronTests.*(neuronTests+1)./((length(ind)-neuronTests-1));
    figure(2),clf;
    plot(neuronTests,  AICtest, '-o', 'Markersize', 7.0, 'LineWidth', 1.5);
    ylabel('AIC^C');
    xlabel('No of hidden neurons');
    box on;
    print('-deps', '-r300', 'exmaple8_3.eps');  
end

% 3 is obtained after executing the above testNeuron loop
numHiddenNeurons = 5;
net = newff(inputvec,targetvec, numHiddenNeurons, {'tansig'});
%Do ANN ensemble if necessary, which is not done in this example
%To do ANN ensemble, use net=init(net) to re-initiate net each time

net.divideFcn='divideblock';
net.divideParam.trainRatio = 0.65;  % Adjust as desired
net.divideParam.valRatio = 0.15;  % Adjust as desired
net.divideParam.testRatio = 0.20;  % Adjust as desired
net.trainFcn = 'trainlm';

Nens=1;
ensemble=cell(Nens, 1);
P = inputvec;
T = targetvec;
% Set up an ensemble of ANN by initializing with random initial weights 
for iens=1:Nens
    net = init(net);
    % Train and Apply Network
    [net,tr] = train(net,inputvec,targetvec);
    trainInd = tr.trainInd;
    valInd = tr.valInd;
    testInd = tr.testInd;
    % Simulate for testing period
    [Y3,~,~,~,info.test.performance] = sim(net,P(:,testInd));
    %simulate for training period
    [Y1,~,~,~,info.train.performance] = sim(net,P(:,trainInd));
    %simulate for validation period
    [Y2,~,~,~,info.validation.performance] = sim(net,P(:,valInd));
    % Get R
    [~,~,info.train.regression] = postreg(Y1,T(:,trainInd),'hide');
    [~,~,info.validation.regression] = postreg(Y2,T(:,valInd),'hide');
    [~,~,info.test.regression] = postreg(Y3,T(:,testInd),'hide');
    disp(['No ' num2str(iens) ' Test error=' num2str(info.test.regression)]);
    ensemble{iens}=net;
end
ssTrain = mapminmax('reverse', Y1, SS);
ssVal = mapminmax('reverse', Y2, SS);
%% Plotting
xi = 1:nrec;
obsval=transpose(ss0(offset+testInd));
ssTest = mapminmax('reverse', Y3, SS);
ssTest = (ssTest.*lambda+1).^(1/lambda);

figure(1), clf;
hold on;
plot(obsval, 'o');
plot(ssTest, '-', 'Color', [0.7, 0.7 0.7], 'LineWidth', 1);
dateaxis('X', 2, alldates(testInd(1)+offset));
hold off;
box on;
xlabel('Time (day)');
ylabel('Q_{transformed}');
legend('Observed', 'ANN', 'Location', 'SouthWest');
legend boxoff;
%% Calculate statistics
E = 1 - sum((ssTest-obsval).^2)/sum((obsval-mean(obsval)).^2)
n = length(ssTest);
MSE =  sum((obsval-ssTest).^2)./n;
K = numHiddenNeurons;
AIC = n*log(MSE)+2*K+2*K*(K+1)/(n-K-1)
RMSE = sqrt(MSE)
basemodelMSE=MSE;

%% sensitivity analysis
doSensitivity=true;
if (doSensitivity)
    inputRange=cell(3,1);
    [DIM1, ~] = size(inputvec);
    inputRange{1} = 3:DIM1; %no precip
    inputRange{2} = [1:3 6:DIM1]; %no tmax
    inputRange{3} = [1:5 8:DIM1]; %no tmin
    for itest=1:3
        Xin = inputvec(inputRange{1},:);
        net = newff(Xin,targetvec, numHiddenNeurons, {'tansig'});
        net.divideFcn='divideblock';
        net.divideParam.trainRatio = 0.65;  % Adjust as desired
        net.divideParam.valRatio = 0.15;  % Adjust as desired
        net.divideParam.testRatio = 0.20;  % Adjust as desired
        net.trainFcn = 'trainlm';


        % Train and Apply Network
        [net,tr] = train(net,Xin,targetvec);
        trainInd = tr.trainInd;
        valInd = tr.valInd;
        testInd = tr.testInd;
        % Simulate for testing period
        [Y3,~,~,~,info.test.performance] = sim(net,Xin(:,testInd));

        qS = mapminmax('reverse', Y3, SS);
        qS = (qS.*lambda+1).^(1/lambda);

        MSE1 =  sum((obsval-qS).^2)./n;
        disp(['Sens Model' num2str(itest) ' ' num2str(MSE1/basemodelMSE)]);
    end
end
    