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
%Dependencies: libsvm-3.11 (included)
%Author: Alex Sun
%Date: 02042013
%Example 8.7
%Simple SVM benchmarking
%==========================================================================
clear all;
addpath(genpath('./libsvm-3.11'));
%Generate training data
incr=0.02;
xleft =-10:incr:-incr;
xright=incr:incr:10;
x=[xleft 0 xright];
y0=[sin(xleft)./xleft 1 sin(xright)./xright];
%reset random seed
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',122361));

%Add noise
y=y0+0.1*randn(size(y0));
%Randomly pick training samples
indx=randperm(length(y)+1);
yt=y0(indx(1:50));
ys=y(indx(1:50));
xs=x(indx(1:50));
%
rngC = mean(ys)-3*std(ys):0.5:mean(ys)+3*std(ys)

figure(1), clf;
plot(x, y0, '-', xs, ys, '+');



%% SVM (LibSVM instruction)
% -s svm_type : set type of SVM (default 0)
% 	0 -- C-SVC
% 	1 -- nu-SVC
% 	2 -- one-class SVM
% 	3 -- epsilon-SVR
% 	4 -- nu-SVR
% -t kernel_type : set type of kernel function (default 2)
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
% 	4 -- precomputed kernel (kernel values in training_set_file)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
% -v n: n-fold cross validation mode
% -q : quiet mode (no outputs)
param.s = 3; 					% epsilon SVR
param.t = 2; 					% RBF kernel
param.nfold = 5;				% 5-fold CV
my = mean(ys);
stdy = std(ys);
param.cset = ...
    1.0:0.1:2.0;   % C
param.gset = 0.1:0.1:1.0;       % range of the gamma parameter
param.p    = 0.2;               % epsilon 
param.e    = 1e-4;              % tolerance

% build model on Learning data
trainvec=xs';
targetvec=ys';

crossval=true;
if (crossval)
bestcv=0;
for j = 1:length(param.gset)
    param.g = param.gset(j);
    for k = 1:length(param.cset)
        param.c=param.cset(k);
        param.libsvm = ['-s ', num2str(param.s), ' -t ', num2str(param.t), ...
            ' -c ', num2str(param.c), ' -g ', num2str(param.g), ...
            ' -p ', num2str(param.p), ' -v ', num2str(5), ' -h 0'];
        cv  = svmtrain(targetvec,trainvec, param.libsvm);
        if (cv >= bestcv)
              bestcv = cv; 
              bestc = param.c; 
              bestg = param.g;
        end
    end
end
    save('example8_6.mat', 'bestcv', 'bestc', 'bestg');
else
    load('example8_6.mat');
end

bestc = 1.5
bestg = 0.12;
param.libsvm = ['-s ', num2str(param.s), ' -t ', num2str(param.t), ...
    ' -c ', num2str(bestc), ' -g ', num2str(bestg), ...
    ' -p ', num2str(param.p)];
model  = svmtrain(targetvec,trainvec, param.libsvm);

testvec = zeros(size(y0'));

inputvec=x';
[ssTest, Acc, projection] = svmpredict(testvec, inputvec, model);

bestc2 = 1.5
bestg2 = 0.12;
param.p    = 0.1;               % epsilon 

param.libsvm = ['-s ', num2str(param.s), ' -t ', num2str(param.t), ...
    ' -c ', num2str(bestc2), ' -g ', num2str(bestg2), ...
    ' -p ', num2str(param.p)];
model2  = svmtrain(targetvec,trainvec, param.libsvm);

inputvec=x';
[ssTest2, Acc, projection] = svmpredict(testvec, inputvec, model2);
%Get the SVs
SVs = full(model.SVs);
indx=[];
for i=1:length(SVs)
    indx = [indx find(xs==SVs(i))];
end
%Figure 8.6 in the Chapter 8
figure(2), clf;
gcolor = [0.165, 0.58, 0.275];
subplot(1,2,1);

hold on;
plot(x, y0, '-', 'LineWidth', 1.5);
plot(x, ssTest, ':', 'LineWidth', 1.5)
plot(xs, ys, '+', 'Color', gcolor, 'LineWidth', 1.5);
plot(xs(indx), ys(indx), 'ro');
hold off;
box on;
xlabel('x');
ylabel('y');
legend('sinc(x)', 'Train Data', 'SVM', 'SV');

%Get the SVs
SVs = full(model2.SVs);
indx=[];
for i=1:length(SVs)
    indx = [indx find(xs==SVs(i))];
end
subplot(1,2,2);
hold on;
plot(xs, ys, '+', 'Color', gcolor, 'LineWidth', 1.5);
plot(x, y0, '-', 'LineWidth', 1.5);
plot(x, ssTest2, ':', 'LineWidth', 1.5);
plot(xs(indx), ys(indx), 'ro');

hold off;
box on;
xlabel('x');
ylabel('y');
%Error
mserr = sqrt(mean((ssTest-y0').^2))
mserr2 = sqrt(mean((ssTest2-y0').^2))