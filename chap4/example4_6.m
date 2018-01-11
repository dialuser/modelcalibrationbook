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

%Author: Alex Sun
%Date: $20111218$
%Example 4.6
%Purpose: Demonstration of Gibbs sampler
%====================================================================
function example4_6()
% prior parameters
m = 9;
s = 3;
a = 3;
b = 0.7;
x = 0:0.02:10;
% note b is inverse of the scale parameter
figure(1), clf;
y=pdf('gam', x, a, 1/b);
plot(x,y);

% set number of samples for burnin and Gibbs
niter = 5000;
nburnin = 500;
ntotal = niter+nburnin;
%true population stat
mu_true = 11;
s_true = 2;
%reset random generator seed
rng = RandStream.create('mt19937ar','seed',2017);
RandStream.setGlobalStream(rng);

n=50; %number of samples
xvec = mu_true+s_true*randn(n,1);
if (exist('example46data.mat', 'file'))
    %use saved results so that experiment can be repeatable
    load('example46data.mat');
else
    save ('example46data.mat', 'xvec');
end
xmean = mean(xvec)
xtot = sum(xvec);
allsamples = zeros(ntotal,2);
% Initial guess for sigma
sigma = 1.5;
for i=1:ntotal
    s12 = 1/(1/s^2+n/sigma^2);
    m1 = (m/s^2+xtot/sigma^2).*s12;    
    %sample mu from N(m1,s1)
    mu = normrnd(m1,sqrt(s12));
    a1 = a+n/2;
    b1 = 0.5*sum((xvec-mu).^2)+b;
    %sample 1/sigma^2 from Gamma(a1,b1)
    lambda = gamrnd(a1,1/b1);
    sigma = sqrt(1/lambda);
    allsamples(i,2) = sigma;
    allsamples(i,1) = mu;
end
%% Do some plotting
allsamples=allsamples(nburnin+1:ntotal,:);
mus = allsamples(:,1);
sigmas = allsamples(:,2);
figure(2), clf;
subplot(2,2,1);
rhist(mus,20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8]);
mean(mus)
var(mus)
xlabel('\mu');
subplot(2,2,2);
rhist(sigmas,20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8]);
mean(sigmas)
var(sigmas)
xlabel('\sigma');
%do trace plot
subplot(2,2,3); plot(mus); ylabel('\mu'); xlabel('Sample No.');
xlim([0 5000]);
subplot(2,2,4); plot(sigmas); ylabel('\sigma'); xlabel('Sample No.');
xlim([0 5000]);

