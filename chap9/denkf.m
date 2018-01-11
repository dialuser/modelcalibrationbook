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

function [ens_mean0, ens_mean, rmse_assi, indx] = ...
    denkf(nrz, nsamples, nparams, truekfield, allsamples)
%% denkf.m
% Implements the deterministic EnKF 
% Reference: Pavel Sakov, Tellus (2008), 
% A deterministic formulation of the ensemble Kalman filter: an 
% alternative to ensemble square root filters
% Last update:20080328
% Fixed bug: use updated head for the next run

global ts;
global errsigma;
global trueindx;
global Tmat;
global h0;

nt = length(ts);
truekfield = log(truekfield);
nstates = 2*nparams;
Amat = zeros(nstates, nrz);
%set initial conditions
Amat(1:nparams,:) = log(Tmat(:,1:nrz));
Amat(nparams+1:nstates, :) = h0(:, ones(nrz,1));
cov_d = errsigma^2.*eye(nsamples);  
rmse_assi = zeros(nt+1, 1);
indx=1:nrz;
%
%assimilation loop starts
%
for it = 1:nt
    dt = ts(it);
    asample = allsamples(:,it);
    for k=1:nrz
        %get initial conditon
        h = Amat(nparams+1:nstates, k);
        %run forward model
        if (it==1)
            T = Tmat(:, k);
        else
            T = exp(Amat(1:nparams, k));
        end
        h = hsolver(h, dt, T);
        Amat(nparams+1:nstates, k) = h;
    end
    Amean = mean(Amat,2);
    if (it==1)
        ens_mean0 = Amean(1:nparams);
        rmse_assi(it) = sqrt(mean((ens_mean0-truekfield).^2))
    end
    %Enkf
    A1 = Amat - Amean(:, ones(nrz,1));
    HA1 = A1(nparams+trueindx, :);
    D1 =  asample - Amean(nparams+trueindx);
    X2 =  inv(HA1*HA1' + (nrz-1)*cov_d);    
    X3 = A1*HA1'*X2; %this is Kalman gain
    %Update mean
    Amean = Amean + X3*D1;
    %generate new variates
    A1 = A1 - 0.5*X3*HA1;
    Amat = Amean(:, ones(nrz,1)) + A1;
    ens_mean = Amean(1:nparams);    
    rmse_assi(it+1) = sqrt(mean((ens_mean-truekfield).^2))    
end