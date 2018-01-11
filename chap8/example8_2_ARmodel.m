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
%Dependencies: Matlab's Financial toolbox (boxcox, dateaxis) and 
%Econometrics toolbox (autocorr function)

%Author: Alex Sun
%Date: 02042013
%Example 8.2
%Purpose: Demonstrate AR process for streamflow forecasting (daily)
%Data: 
%Streamflow: USGS gauge 08033500, Neches Rv nr Rockland, TX
%http://waterdata.usgs.gov/nwis/nwisman/?site_no=08033500
%3,636 square miles
%Duration 1948/1/1-2003/12/31
%Other data: MOPEX 
%ftp://hydrology.nws.noaa.gov/pub/gcip/mopex/US_Data/Us_438_Daily/08055500.
%dly

%Results:
% Usage: ORDERS = STRUC(NA_RANGE,NB_RANGE,NK_RANGE)
% Discrete-time IDPOLY model: A(q)y(t) = e(t)                           
% A(q) = 1 - 1.457 q^-1 + 0.574 q^-2 - 0.1625 q^-3 + 0.05154 q^-4       
%         - 0.004798 q^-5 + 0.008775 q^-6 + 0.01038 q^-7 - 0.00251 q^-8 
%                                        - 0.005677 q^-9 - 0.01108 q^-10
% Estimated using ARMAX from data set z                                 
% Loss function 0.000329073 and FPE 0.000329474                         
% Sampling interval: 1      
%==========================================================================
%% Load daily streamflow, Q in [cfs]
clear all
fid = fopen('08033500_Q.dat', 'r');
res = textscan(fid, '%s %f');
fclose(fid);
alldates=res{1}(:); %date strings
ss0 = res{2}(:);
% Transform the streamflow data using box-cox transform
[ss1, lambda] = boxcox(ss0);
ss = ss1;

%% Normalize data
nrec = length(ss);
%scaling
XMIN = 0.01
XMAX = 0.99;

[ss, SS] = mapminmax(ss', XMIN, XMAX);    

%% Construct AR model
data=iddata(ss',[], 1);
%This corresponds to end of 1992
nTrain = 16437;
%data used for training
ze = data(1:nTrain);
%data used for validation
zv = data(nTrain+1:end);
%Get the best orders
V = arxstruc(ze,zv,struc(1:10));
nn = selstruc(V,0);

%order of the autogressive variable
na =nn(1);
%order of the moving averaging 
nc = 0;
ssmodel=armax(ze, 'na', na, 'nc', nc)

[yh] = compare(zv, ssmodel,1);
ssTest =yh.OutputData; 

%for k-step ahead prediction (if performKStep is set to true)
performKStep=false;
if (performKStep)
    %Perform n-step-ahead prediction
    K=6;
    yh = predict(ssmodel, zv, K);
    %plot
    predict(ssmodel, zv, K);
    compare(zv, ssmodel, K)
    ssTest =yh.OutputData;
end
% one-step-ahead prediction without calling matlab function
forecasting=false;
if (forecasting)
    %do AR(p) for future values
    %get AR(p) coefficients from a trainined model
    paramAR = -ssmodel.a(2:end);
    yn = ze(end-24:end).OutputData;
    K=10;
    ssTest3 = zeros(1,K);
    for k=1:K
        yn1=paramAR*yn;
        ssTest3(k)=(yn1);

        yn(1:end-1)=yn(2:end);
        yn(end)=yn1;
    end
    origdata = ss(nTrain+1:nTrain+K);
    scatter(ssTest3, origdata )
    xlim([0 1]);ylim([0 1]);
end

%% Plotting
close all;
figure(1), clf;
autocorr(ss, 365);

figure(2);
% scatter plot on transformed data
%Transform data back in the original space
ssTest = mapminmax('reverse', ssTest, SS);
obsval = ss0(nTrain+1:end);
ssTest = (ssTest.*lambda+1).^(1/lambda);

subplot(1,2,1);
scatter(log(ssTest), log(obsval));
box on;
xlabel('Predicted, log(Q)');
ylabel('Observed, log(Q)');
subplot(1,2,2);
qqplot(ssTest, obsval);
set(gca, 'YTick', 0:1e4:5e4);
xlabel('Predicted, Q');
ylabel('Observed,  Q');
box on;
%
figure(3);
hold on;
xi = 1:length(obsval);
plot(xi, obsval, 'o');
plot(xi, ssTest, 'r', 'LineWidth', 1.0, 'Color', [0.7,0.7,0.7]);
xlabel('Date ');
ylabel('Q (cfs)');
box on;
xlim([0 length(xi)])
dateaxis('X', 12, alldates(nTrain+1));
legend('observed', 'AR model', 'Location', 'SouthWest');

%Calculate NSE
E = 1 - sum((ssTest-obsval).^2)/sum((obsval-mean(obsval)).^2)
MSE = mean(sum(ssTest-obsval).^2)

