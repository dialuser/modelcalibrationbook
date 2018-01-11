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
%Example 2.6: linearization
%Author: Alex Sun
%Purpose: Linearization of kinematic wave equation
%Data: example 9.4 in Chow W.T. et al. 1988
%====================================================================

function []=example2_6()
%Boundary condition
hb=[2000 2000 2000 2000 2000 2250 2500 2750 3000 3250 3500 3750 ...
4000 4250 4500 4750 5000 5250 5500 5750 6000 5750 5500 5250 5000 ...
4750 4500 4250 4000 3750 3500 3250 3000 2750 2500 2250 2000 2000 ...
2000 2000 2000];
%Discretization
dx=1000;
x=dx:dx:24000;
dt=3;
t=0:dt:120;
dt=3*60; %convert to seconds
N=length(x);
h0=zeros(size(x))+2000;
q=0;
beta=0.6;
n=0.035;
P=200;
Sf=0.01;
sfRange=0.005:0.001:0.015;
%Sensitivity to energy slope, Sf
res = zeros(size(sfRange));
for k=1:length(sfRange);
    alpha=getAlpha(n,P,sfRange(k));
    h=getSolution(t,N,alpha,beta,h0,hb,q,dt,dx);
    res(k) = h(10,3);
end
figure(1), clf;
subplot(1,2,1);
plot(t,hb, 'r', 'LineWidth', 1.5);
xlabel('t (min)');
ylabel('Discharge (ft^3/s)');
hold on;
plot(t,h(:,2), '--', 'LineWidth', 1.5);
plot(t,h(:,10), '-.', 'LineWidth', 1.5);
legend('Inflow', 'x=2\Delta', 'x=10\Delta');
legend boxoff;
box on;
hold off;
xlim([0 120]);

subplot(1,2,2);
plot(sfRange, res);
xlabel('S_0');
ylabel('Discharge (ft^3/s)');
%Get true observations
Sf=0.01;
alpha=getAlpha(n,P,Sf);
htrue=getSolution(t,N,alpha,beta,h0,hb,q,dt,dx);
obsloc=12;
trange=10:4:30;
obsRange=[trange' zeros(length(trange),1)+obsloc];
nobs = size(obsRange,1);
trueobs=getObs(htrue,obsRange);
RandStream.setGlobalStream... 
     (RandStream('mt19937ar','seed',19319));
obserr= 1*randn(nobs,1);
bvec=trueobs+obserr;
figure(2), clf;
subplot(1,2,1);
hold on;
plot(t, htrue(:,12), '-', 'LineWidth', 1.5);
plot(trange*3-3, bvec, 'o');
box on;
hold off;
xlabel('t (min)');
ylabel('Dicharge (ft^3/s)');

%Start linearization process
%Get sensitivity matrix
%Initial guess
Sf0=0.005;
epsil = 1e-2; %convergence criterion
for k=1:100
    alpha=getAlpha(n,P,Sf0); 
    h=getSolution(t,N,alpha,beta,h0,hb,q,dt,dx);
    Sf1=Sf0*(1+0.05);
    alpha1=getAlpha(n,P,Sf1);
    h1=getSolution(t,N,alpha1,beta,h0,hb,q,dt,dx);
    Jall=(h1-h)./(Sf1-Sf0);
    hm = getObs(h,obsRange);
    J=getObs(Jall,obsRange);
    if (k==1)
        subplot(1,2,2);
        hold on;
        plot(t, Jall(:, obsloc), 'LineWidth', 1.5);
        plot(trange*3-3, J, 'o', 'LineWidth', 1.5);
        xlabel('t (min)');
        ylabel('J=dQ/dS_0 (ft^3/s)');
        box on;
        hold off;
    end
    b=bvec-hm+J*Sf0;
    theta = J\b
    if (abs((theta-Sf0)/Sf0)<epsil) 
        disp(theta);
        disp(k);
        break;
    end
    Sf0=theta;
end
function h=getSolution(t,N,alpha,beta,h0,hb,q,dt,dx)
%Solve water level at all nodes using the finite-difference method    
h=zeros(length(t), N);
for i=1:length(t)
    for j=1:N
        if (j==1)
            meanh = alpha*beta*(0.5*(h0(j)+hb(i)))^(beta-1);
            h(i,j)=(dt/dx*hb(i)+meanh*h0(j)+q*dt)/(dt/dx+meanh);
        else
            meanh = alpha*beta*(0.5*(h(j-1)+h0(j)))^(beta-1);
            h(i,j)=(dt/dx*h(i,j-1)+meanh*h0(j)+q*dt)/(dt/dx+meanh);            
        end
    end
    h0=h(i,:);
end
function res=getAlpha(n,P,Sf)
    res=(n*P^(2/3)/(1.49*sqrt(Sf)))^(3/5);
function res=getObs(amat, obsRange)
    nobs = size(obsRange,1);
    res=[];
    for kk=1:nobs
        res=[res; amat(obsRange(kk,1), obsRange(kk,2))];
    end
