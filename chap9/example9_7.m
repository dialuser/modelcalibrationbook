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
%MODEL CALIBRATION AND PARAMETER ESTIMATION
%Author: Alex Sun
%Revision $20130625$
%Dependencies:
%   hsolver_ex2.m
%   hsolver_ex3.m
%   ekfukf toolbox 
%        Optimal filtering with Kalman filters and smoothers – a
%        Manual for Matlab toolbox EKF/UKF
%        Jouni Hartikainen and Simo Särkkä
%        Department of Biomedical Engineering and Computational Science,
%        Helsinki University of Technology,
%        jmjharti@lce.hut.fi, simo.sarkka@hut.fi
%==========================================================================
% Example9_7.m, Example9.6 in the book
% Description: Demonstration of UKF
% One-dimensional groundwater flow subject to constant-head BC H0 and 
%   constant-flux boundary  
% In this example, transmissivity is the uncertain parameter to be
% estimated. 
% Note hsolver_ex2 and hsolver_ex3 do essentially the same thing; however
% they provide different functional interface. I could have combined into a
% single function
%
%% Set up model parameters
close all;
addpath ('./ekfukf');

global q S;
global nx dx;
global H0 qL;

H0=10; 
qL=1e-3;
nx=100;

S  = 1e-5;
Ym = -4;   %actual value
Y0 = -4.4; %initial guess

%True transmissivity
T = exp(Ym); 
%Initial conditions
h0 = zeros(nx,1) + 5.0;
%Pumping rates one well at 70-th cell
q = zeros(1,nx); 
q(70) = -2e-4;
%Time steps
dt = 0.05;
tmax = 1.0;
ts = dt:dt:tmax;
nt = length(ts);
dx = 1.0;
nx = 100;
%% Generate synthetic measurements
trueindx = 10:10:nx-10; 
nsamples = length(trueindx);
allsamples = zeros(nsamples, nt);
 
figure(1), clf;
hold on;
h=h0;
for i=1:nt
    [h Amat bvec] = hsolver_ex2(h, dt, log(T));
    allsamples(:,i) = h(trueindx);    
    plot(trueindx, h(trueindx), 'o', 'LineWidth', 1.0);       
end
%% Unscented KF
% Set up Amat
h=h0;
% Set up error matrix
R = 1e-2*eye(nsamples,nsamples);
sigma_eta = 1e-2;
Q_eta = sigma_eta^2*eye(nx,nx);
sigmav= 1e-2;
Q_v = sigmav^2*eye(1,1);
G = zeros(nsamples, nx+1);
for i=1:nsamples
    G(i, trueindx(i)) = 1;
end
P_h=1e-6*eye(nx,nx);

varY = 0.3;
P_Y = varY;
Y=Y0;

disp(['Initial T= ' num2str(T)]);
P = [P_h zeros(nx,1); zeros(1,nx) P_Y];
Q = [Q_eta zeros(nx, 1); zeros(1,nx) Q_v];

fhandle = @hsolver_ex3;
y=[h0;Y];
params={dt};
simiVal=zeros(nt,1);
PyVal = zeros(nt,1);

for i=1:nt
    [y,P] = ukf_predict1(y,P, fhandle, Q, params);
    [y,P] = ukf_update1(y,P,allsamples(:,i),G,R);
    T=exp(y(end))
    simiVal(i) = T;
    h=y(1:nx);   
    PyVal(i) = P(end,end);
    if (i==1)
        plot(h, 'c-.', 'Linewidth', 1.5);
    elseif (i==nt)
        plot(h, 'g-.', 'Linewidth', 1.5);
    else
        plot(h, '-.', 'Linewidth', 1.5, 'Color', [0.8 0.8 0.8]);
    end        
end
hold off;
box on;
xlabel('x, [L]');
ylabel('Head h(t), [L]');

%load ekf.mat for comparison
load('ekf.mat');

figure(2), clf;
hold on;
plot(ts, simiVal, 'x', 'Linewidth', 1.5);
plot(ts, simiValEKF, 'o', 'Linewidth', 1.5);
xlabel('Assimilation Time, [T]');
ylabel('K, [L/T]');
xlim([0 tmax]);
box on;
hold off;
legend('UKF', 'EKF', 'Location', 'best');

figure(3), clf;
hold on;
plot(ts, PyVal, '^', 'Linewidth', 1.5);
plot(ts, PyEKF, 'o', 'Linewidth', 1.5);
xlabel('Assimilation Time, [T]');
ylabel('P_Y');
xlim([0 tmax]);
box on;
hold off;
