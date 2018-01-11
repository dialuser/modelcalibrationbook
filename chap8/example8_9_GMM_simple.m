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
%Code dependencies: GPML (included)

%Author: Alex Sun
%Example 8.9 GP demonstration
%Demonstrate sinc(x) fitting using GPML
%====================================================================
clear all, close all;
addpath(genpath('./gpml-matlab'));

%Generate training data
incr=0.02;
xleft =-10:incr:-incr;
xright=incr:incr:10;
x=[xleft 0 xright]';
y0=[sin(xleft)./xleft 1 sin(xright)./xright]';
%Reset random seed
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',122361));

%Add noise to "true" data
stdY = 0.1;
n = 50;
y=y0+stdY*randn(size(y0));
%Randomly pick training samples
indx=randperm(length(y)+1);
yt=y0(indx(1:n));
ys=yt+stdY*randn(size(yt));
xs=x(indx(1:n));
%% GPML 
% specify prior model for mean, cov, and likelihood
meanfunc = {@meanConst};hyp.mean = 0.0;
%ell=length scale; 
covfunc = {'covSEisoU'};ell=1.0; hyp.cov= log(ell);   % isotropic Gauss no scale
%sn = noise standard deviation!
likfunc = @likGauss; sn = stdY; hyp.lik = log(sn);

K = feval(covfunc{:}, hyp.cov, xs);
mu = feval(meanfunc{:}, hyp.mean, xs);

%confidence for the mean function
hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, xs, ys);
[my s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xs, ys, x);
figure;
%95% confidence interval
f = [my+2*sqrt(s2); flipdim(my-2*sqrt(s2),1)];
fill([x; flipdim(x,1)], f, [7 7 7]/8)
hold on; 
plot(x, y0, '-', 'LineWidth', 1.5);
plot(x, my, '--', 'LineWidth', 1.5);
plot(xs, ys, '+');
xlabel('x');
ylabel('y');
legend('Bounds', 'sinc(x)', 'GP', 'Train data');
hold off;
rmse = sqrt(mean((my-y0).^2));
disp(rmse);


  