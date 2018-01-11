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
%Dependencies:
%   hsolver_ex1.m
%   ekfukf toolbox (Jouni Hartikainen and Simo Särkkä)

%Author: Alex Sun
%Revision $20130625$
%
%==========================================================================
function example9_4
%Example 9.4 and 9.5
%Description: Demonstration of classic KF
%One-dimensional groundwater flow subject to constant-head BC, [H0, H1]. 
%In this example, H0 is uncertain and is assumed as a Gaussian random
%variable and our goal is to identify H0 by using joint state and
%'parameter' KF
%
%Joint state and parameter estimation
%Augmented variable y = [h; H0]
%State equation: h_t = A*h_(t-1) + Bu + eta
%Parameter equation: H0 = H0 + nu
%Measurement equation: z_t = G*h_(t-1)+epsilon


addpath ('./ekfukf');

%% Set up model parameters
global q S;
global nx dx;
global H0 HL;
global a1;

H0=10; %This is the "true" left-hand boundary condition
HL=1;
%discretization along x direction
nx=100;
%specific storage
S = 1e-5;
Ym = -4;
%hydraulic conductivity
T = exp(Ym); 
%Initial conditions
h0 = zeros(nx,1) + 5.0;
%Pumping rates, one well located at the 70-th cell
q = zeros(1,nx); 
q(70) = -2e-4;
%Time discretization
dt = 0.1;
tmax = 0.5;
ts = dt:dt:tmax;
nt = length(ts);
dx = 1.0;
nx = 100;
%% Generate measurements
trueindx = 10:10:nx-10; 
nsamples = length(trueindx);
allsamples = zeros(nsamples, nt);

figure(1), clf;
subplot(1,2,1), hold on;
h=h0;
for i=1:nt
    [h Amat bvec cvec] = hsolver_ex1(h, dt, T);
    truehfield(:,i) = h;
    %plot(h, 'LineWidth', 1.5)
    allsamples(:,i) = h(trueindx);
    plot(trueindx, h(trueindx), 'o', 'LineWidth', 1.0);    
end

%% Start KF
%Initial guess of boundary condition
H_guess = 8;
H0=H_guess;
disp(['Initial guess of boundary condition ' num2str(H0)]);

% reset random generator seed to make results reproducible
rng = RandStream.create('mt19937ar','seed',69999);
RandStream.setGlobalStream(rng);

h=h0;
%define augmented state vector
y=[h; H0];
% Set up measurement error matrix
R = 0.01*eye(nsamples,nsamples);
sigma_eta= 0.1;
% model-state error matrix
Q_eta = sigma_eta^2*eye(nx,nx);
% parameter error 
sigma_v = 0.1;
Q_v = sigma_v^2*eye(1,1);
% set up measurement operator
G = zeros(nsamples, nx+1);
for i=1:nsamples
    G(i, trueindx(i)) = 1;
end
P_h=zeros(nx,nx)+0.0;
P_H0 = 0.1;
P = [P_h zeros(nx,1); zeros(1,nx) P_H0];
[h Amat bvec cvec]=hsolver_ex1(h, dt, T);

Amat = inv(Amat);
Q = [Q_eta zeros(nx, 1); zeros(1,nx) Q_v];
A = [a1.*Amat Amat*cvec; zeros(1,nx) 1];
B = [Amat zeros(nx,1); zeros(1,nx+1)];
u = [bvec; 0];

MM = zeros(nx+1, nt);
PP = zeros(nx+1, nx+1, nt);
% KF iterations
simiVal=zeros(length(nt), 1);
for i=1:nt
    [y,P]= kf_predict(y, P, A, Q, B, u);
    asample = allsamples(:,i);
    [y,P]= kf_update(y, P, asample, G, R); 
    if (i==1)
        plot(y(1:nx), 'c:', 'Linewidth', 1.5);
    elseif (i==nt)
        plot(y(1:nx), '-', 'Linewidth', 1.5, 'Color', [153 76 0]/255);
    else
        plot(y(1:nx), 'g-.', 'Linewidth', 1.5);
    end
    MM(:,i) = y;
    PP(:,:,i) = P;
    H0=y(end)
    simiVal(i) = H0;
end
hold off;
xlabel('x (L)');
ylabel('Head, h(t) (L)');
box on;
subplot(1,2,2);

%add the initial guess for plotting
ts = [0.0 ts];
simiVal =[H_guess simiVal];
plot(ts, simiVal, 'x', 'Linewidth', 1.5);
xlabel('Assimilation Time (T)');
ylabel('H_0 (L)');
xlim([0 tmax]);
ylim([H_guess  10]);
% KF smoother
[SM,SP] = rts_smooth(MM,PP,A,Q);