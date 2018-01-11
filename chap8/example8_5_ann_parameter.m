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
%Dependencies: Matlab ANN toolbox, sampling(included), global optimization
%Example 8_5.m
%Demonstrate function approximation using ann
%====================================================================
close all;
myfunc=@(x, p) p(1)*sin(p(2)*x)+p(3)*sin(p(4)*x);
%Generate true value
ptrue=[0.25 -0.51 0.23 1.1]';
x=0:0.1:2*pi();
y_true=myfunc(x,ptrue);
%define parameter intervals
pint=[0.1 0.5
      -1  0
      0   0.5
      0.5 1.5];
%% Do ann approximation
numNeurons = 15;
trainnet=true;
if (trainnet)
    %Generate training samples
    addpath('./sampling');
    RandStream.setGlobalStream ... 
         (RandStream('mt19937ar','seed',123375));
    % define random dimension
    N = size(pint, 1);
    nrz=200;
    allSamples = lhsu(zeros(1,N), zeros(1,N)+1, nrz);    
    inputvec=zeros(N, nrz);
    targetvec=zeros(length(x), nrz);
    for irz=1:nrz
        rnd = allSamples(irz,:);
        tempvec = pint(:,1)+(pint(:,2)-pint(:,1)).*rnd';
        inputvec(:,irz) = tempvec;
        targetvec(:,irz) = myfunc(x, tempvec);
    end

    net = newff(inputvec,targetvec, numNeurons);
    %Do ANN ensemble if necessary, which is not done in this example
    %To do ANN ensemble, use net=init(net) to re-initiate net each time

    net.divideFcn='divideblock';
    net.divideParam.trainRatio = 0.65;  % Adjust as desired
    net.divideParam.valRatio = 0.15;  % Adjust as desired
    net.divideParam.testRatio = 0.20;  % Adjust as desired
    net.trainFcn = 'trainlm';
    % Train and Apply Network
    [net,tr] = train(net,inputvec,targetvec);
    save('example8_5.mat', 'net');
else
    load('example8_5.mat');
end
%% forward model comparison with true model
inpoint = ptrue;
outpoint = sim(net,inpoint);
%% inversion using GA
%setup ga solver
options = gaoptimset('Generations',100);
strm = RandStream('mt19937ar','Seed',16525);
RandStream.setGlobalStream(strm);
options.PopulationSize = 30;
options.ParetoFraction = 0.35;
options.MigrationFraction=0.25;
options.TolFun=1e-5;
tstart=tic;
annhandle = @(pvec) sqrt(mean((y_true'-sim(net,pvec')).^2));
[gaSol,fval,exitflag,output, final_pop] = ga(annhandle, 4, [],[],[],[], ...
                                        pint(:,1), pint(:,2), [], options);
telapsed = toc(tstart)                                    
disp(gaSol);  
figure(1), clf;
plot(x, y_true, 'o', x, outpoint, 'r-', x, myfunc(x, gaSol), 'k-.', 'LineWidth', 1.5);
xlabel('x');
ylabel('f(x)');
box on;
legend('True', 'ANN Approximation', 'ANN-GA', 'Location', 'Best');
legend boxoff;
xlim([0 2*pi()])
MSE = sqrt(mean(sum(outpoint-y_true').^2))
