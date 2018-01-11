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
%Revision $20130630$
%Dependencies:
%   hsolver.m
%Example9_8.m
%Purpose: One-d test problem
%
% nrz, number of realizations in the ensemble
% nx, dimension of discretization along the one-d column
% dx, size of each cell
% hR, hL, dirichlet boundary values
% h0, initial guess of state vector
% Q, collection of all sink/source terms
% errsigma, measurement error standard deviation
% Ym, sigmaY, iscaleY: mean, stddev, and correlation len of log(conductivity)
% Tmat, matrix random conductivity fields in columns
% allsamples, collection of measurements for all time steps
% S, specific storage [1/L]
%To use the code:
% First, generate results for both EnKF and DEnKF by setting the rerun
% variable to true and run it once for EnKF and DEnKF by setting value to
% smethod.
% Then, set rerun to false and run the code to plot results
%==========================================================================
%% Gloal Parameters
clear all; close all;
global Q S;
global nx dx;
global ts;
global errsigma;
global trueindx;
global Tmat;
global h0;
rerun = false; %set false to use previously saved results, set true to regen

%% Problem Setup 
nrz = 500;
% reset random generator seed to make results reproducible
rng = RandStream.create('mt19937ar','seed',99);
RandStream.setGlobalStream(rng);

hL = 10;
hR = 1;
dx = 1.0;
nx = 100;
xloc = 0:dx:(nx-1)*dx;
%Initial and boundary conditions
h0 = zeros(nx,1);
h0(1) = hL; h0(nx) = hR;
dh = (hL-hR)/((nx-1)*dx);
h0(2:nx-1) = hL-dh:-dh:hR+dh;
iplot = 0;

Q = zeros(1,nx);
Q(30) = -3e-2;
Q(60) = 3e-2;
S = 1e-5; %specific storage
% Define Assimilation Times
dt = 0.1;
tmax = 0.6;
ts = zeros(ceil(tmax/dt),1)+dt;
nt = length(ts);

%% Parameters for generating single-modal gaussian conductivity fields
%mean log(conductivity)
Ym = -4;
%deviation of log(conductivity)
sigmaY = 0.5;
%integral scale
iscale = 5.0; 
Tmat = kgen(Ym, sigmaY, nx, xloc, iscale, nrz+1);

%%
% Set the true field
truekfield = Tmat(:,1);
% Set the nrz ensemble 
Tmat = Tmat(:,2:nrz+1);
%Define observation locations
%trueindx = 5:5:nx-5; 
trueindx = 10:10:nx-10; 
%trueindx = 2:2:nx-2; 
nsamples = length(trueindx);
allsamples = zeros(nsamples, nt);
errsigma = 1.0e-2;
% Solve for truefield
figure(1), clf;
title('Evolution of True State');
h = h0;
truehfield = zeros(size(h,1), nt);
hold on;
for i=1:nt
    h = hsolver(h, dt, truekfield);
    truehfield(:,i) = h;
    if i==1
        plot(h, '--', 'LineWidth', 1.5);
    else
        plot(h, 'LineWidth', 1.5);
    end
    allsamples(:,i) = h(trueindx);    
end
hold off;
xlabel('Distance x, [L]');
ylabel('Head h, [L]');
box on;
save 'trueh.mat' truehfield;
%% Assimilation
nparams = nx;
if (rerun)
    smethod = 'denkf';
    if (strcmp(smethod, 'denkf'))
      [ens_mean0, ens_mean, rmse_assi, indx] = ...
        denkf(nrz, nsamples, nparams, truekfield, allsamples);    
      denkf_mean0 = ens_mean0;
      denkf_mean = ens_mean;
      denkf_rmse = rmse_assi;
      save 'denkf.mat' 'denkf_mean0' 'denkf_mean' 'denkf_rmse'; 
    elseif (strcmp(smethod, 'enkf'))
      [ens_mean0, ens_mean, rmse_assi, indx] = ...
        enkf(nrz, nsamples, nparams, truekfield, allsamples);    
      enkf_mean0 = ens_mean0;
      enkf_mean = ens_mean;
      enkf_rmse = rmse_assi;
      save 'enkf.mat' 'enkf_mean0' 'enkf_mean' 'enkf_rmse'; 
    end
    save 'truekfield.mat' xloc truekfield
    return
else
    load('enkf.mat');
    load('denkf.mat');
end

%% Display results
figure;
subplot(1,2,1);
hold on;
plot(xloc, log(truekfield), '-', xloc, enkf_mean0, '--', ...
    xloc, enkf_mean, 'ro', 'LineWidth', 1.5);
plot(xloc, denkf_mean, 'g^', 'LineWidth', 1.5);
xlabel('Distance x, [L]');
ylabel('Y=lnK');
legend('True', 'Initial Mean', 'EnKF', 'DEnKF', 'Location', 'best');
hold off;
box on;
subplot(1,2,2);
hold on;
plot(0:dt:tmax, enkf_rmse, 'ro-', 'LineWidth', 1.5);
plot(0:dt:tmax, denkf_rmse, 'g^-', 'LineWidth', 1.5);
xlabel('Time [T]');
ylabel('RMSE');
legend('EnKF', 'DEnKF');
xlim([0 tmax]);    
hold off;
box on;