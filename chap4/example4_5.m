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

%Example 4.5 
%Metropolis Hastings sampling
%====================================================================
function example4_5()
% sigmaq, std dev of the Gaussian proposal
global sigmaq;
N=2000;
%reset random generator seed
rng = RandStream.create('mcg16807', 'Seed',0);
RandStream.setGlobalStream(rng);
%Effect of sigmaq
figure(1), clf;
h=subplot(1,3,1);
[samples1,accrate1, funvals] = runMH(5,10,N);
h=subplot(1,3,2);
[samples2,accrate2, funvals] = runMH(5, 0.25,N);
h=subplot(1,3,3);
[samples3,accrate3, funvals] = runMH(5, 80,N);

figure(2), clf;
%Effect of N
h=subplot(1,2,1);
[samples4,accrate4, funvals] = runMH(5, 10,1000);
h=subplot(1,2,2);
[samples5,accrate5, funvals] = runMH(5, 10,20000);
%Effect of x0
figure(3), clf;
h=subplot(1,2,1);
[samples4,accrate4, funvals] = runMH(5, 10,N);
h=subplot(1,2,2);
[samples5,accrate5, funvals] = runMH(10, 10,N);

%Plot trace
figure(4), clf;
subplot(1,3,1);
plot(samples1);
xlabel('Sample No.');
ylabel('\theta');
subplot(1,3,2);
plot(samples2);
xlabel('Sample No.');
subplot(1,3,3);
plot(samples3)
xlabel('Sample No.');
function [allsamples,accptrate, actualfunc]=runMH(x0, sigmaq,N)
global testfun;
testfun = 2;
if (testfun==1)
    x0=5;
    sigmaq=8;
    xrng=-20:0.5:30;
    lb=-20;
    ub=30;
    Nbins=18;
elseif (testfun==2)
    xrng=-10:0.5:30;
    lb=-10;
    ub=30;    
    Nbins=25;
end
Nburnin=100;
allsamples=zeros(N,1);
naccept=0;
xi=x0;
%burn in
for i=1:Nburnin
    [xi, iaccept] = MHSample(xi, sigmaq);
end
%MH sampling
for i=1:N
    [xi, iaccept] = MHSample(xi, sigmaq);
    naccept=naccept+iaccept;
    allsamples(i)=xi;
end
%Plot statistics
accptrate= naccept/N;
c = quadl(@p, lb,ub);
hold on;
rhist(allsamples, Nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8,0.8,0.8],'EdgeColor','w');
actualfunc=p(xrng)./c;
plot(xrng, actualfunc, 'r', 'LineWidth', 2);
hold off;
box on;
ylabel('PDF');
xlabel('\theta');
%%
function [xnew, iaccept]=MHSample(xi, sigmaq)
%Input variables,
%x, current sample
%sigmaq, standard deviation of the proposal
    
%sample uniform distribution
u = rand();
%sample candidate from normal pdf, N(xi, sigmaq)
xc = xi+randn()*sigmaq;
%calculate ratio
r = p(xc)*q(xi, xc, sigmaq)/(p(xi)*q(xc,xi,sigmaq));
alpha=min([r,1]);
if (u<=alpha)
    xnew = xc;
    iaccept=1;
else
    xnew = xi;
    iaccept=0;
end
%% proposal function
function res=q(xc,x, sigmaq)
res=normpdf(xc, x, sigmaq);

%% target function
function res=p(x)
global testfun;
if (testfun==1)
    %a 'common' bimodal
    res=0.3*exp(-0.2*x.^2)+0.7*exp(-0.2*(x-10).^2);
elseif(testfun==2)
    %a 'distorted' trimodal
    %res= exp(-x.^3).*(2+sin(x)+sin(2*x));
    res=(exp(-(x-5).^2./4)+0.3*exp(-(x-15).^2/10)...
    +0.5*exp(-(x-10).^2./18));
end