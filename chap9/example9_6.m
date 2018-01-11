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
%Revision $20130627$
%Dependencies:
%   hsolver_ex2.m
%   ekfukf toolbox by
%        Optimal filtering with Kalman filters and smoothers – a
%        Manual for Matlab toolbox EKF/UKF
%        Jouni Hartikainen and Simo Särkkä
%        Department of Biomedical Engineering and Computational Science,
%        Helsinki University of Technology,
%        jmjharti@lce.hut.fi, simo.sarkka@hut.fi
%==========================================================================
% Example 9.6 <flux version>
% Description: Demonstration of EKF
% One-dimensional groundwater flow subject to constant-head BC, [H0, H1]. 
% In this example, transmissivity is the uncertain parameter to be
% estimated. 
% Joint estimation of head and transmissivity using EKF
%
% This is a version that one-end is flux boundary
% The assimilation rate is relatively fast compared to the constant-head
% case

%% Set up model parameters
function example9_6
close all;
addpath ('./ekfukf');

global q S;
global nx dx;
global H0 qL;
global a1; 

H0=10; 
qL= 1e-3;
nx= 100;
S = 1e-5;
Ym = -4;
Y0 = -4.4; %initial guess
%Transmissivity
T = exp(Ym);
%Initial condition
h0 = zeros(nx,1) + 5.0;
%Pumping rates one well at 70-th cell
q = zeros(1,nx); 
q(70) = -2e-4;
%Time steps
icase = 1;
switch icase
    case 1
        %high-frequency DA
        dt = 0.05;
        tmax = 1.0;
    case 2
        %low-frequency DA
        dt = 0.1;
        tmax = 1.0;
end
ts = dt:dt:tmax;
nt = length(ts);
dx = 1.0;
nx = 100;
%% Generate synthetic measurements
trueindx = 10:10:nx-10; 
nsamples = length(trueindx);
allsamples = zeros(nsamples, nt);
 
figure(1), clf;
subplot(1,2,1); 
hold on;
h=h0;

for i=1:nt
    [h Amat bvec] = hsolver_ex2(h, dt, log(T));
    allsamples(:,i) = h(trueindx);    
    plot(trueindx, h(trueindx), 'o', 'LineWidth', 1.0);       
end
[simiVal, PyVal] = doEKF

subplot(1,2,2);
plot(ts, simiVal, 'x', 'Linewidth', 1.5);
xlabel('Assimilation Time, [T]');
ylabel('K, [L/T]');
xlim([0 tmax]);

figure(2), clf;
plot(ts, PyVal, '^', 'Linewidth', 1.5);
xlabel('Assimilation Time, [T]');
ylabel('P_Y');
xlim([0 tmax]);



    function [simiVal, PyVal] = doEKF
    %% Extended KF
    % Set up Amat
    h=h0;
    % Set up error matrices
    R = 1e-2*eye(nsamples,nsamples);
    sigma_eta = 1e-2;
    Q_eta = sigma_eta^2*eye(nx,nx);
    sigmav= 0.3;
    Q_v = sigmav^2*eye(1,1);
    G = zeros(nsamples, nx+1);
    for ii=1:nsamples
        G(ii, trueindx(ii)) = 1;
    end
    P_h  =1e-6*eye(nx,nx);
    varY = 0.3;
    P_Y= varY;
    Y = Y0;
    disp(['initial T=', num2str(exp(Y))]);
    P = [P_h   zeros(nx, 1); zeros(1,nx) P_Y];
    Q = [Q_eta zeros(nx, 1); zeros(1,nx) Q_v];
    simiVal = zeros(nt,1);
    PyVal = zeros(nt,1);
    for ii=1:nt
        [ho Amat bvec A1]=hsolver_ex2(h, dt, Y);
        Amat = inv(Amat);
        %derivative w/ T
        dhdT = a1*Amat*A1*Amat*ho;
        A = [a1*Amat dhdT; zeros(1,nx) 1];
        y0=[ho; Y]; %predicted state
        y =[h;  Y]; %previous state
        [y,P] = ekf_predict1(y, P, A, Q, y0);
        [y,P] = ekf_update1(y,P,allsamples(:,ii),G,R);
        Y=y(end);
        simiVal(ii) = exp(Y);
        PyVal(ii) = P(end,end);
        disp([y(end) P(end,end)])
        h=y(1:nx);
        if (ii==1)
            plot(h, 'c-.', 'Linewidth', 1.5);
        elseif (ii==nt)
            plot(h, 'g-.', 'Linewidth', 1.5);
        else
            plot(h, '-.', 'Linewidth', 1.5, 'Color', [0.8 0.8 0.8]);
        end    
    end

    box on;
    xlabel('x, [L]');
    ylabel('Head h(t), [L]');
    hold off;
    
    %save results for comparison with UKF
    simiValEKF = simiVal;
    PyEKF = PyVal;
    save('ekf.mat', 'simiValEKF', 'PyEKF');
    end
end
