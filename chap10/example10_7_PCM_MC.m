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
%Example 10.7
%This is needed by example10_7.m
%Purpose: driver for testing the one-dimensional problem by using MC
%set itype==1 to generate mc_normal.mat
%set itype==2 to generate mc_uniform.mat
%==========================================================================
%Parameters
clear all;
global Q S;
global nx dx;
global ts;
global h0;
%
nrz =1;
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',123374));

nx = 100;
hL = 100;
thick = 10;
hR = 90;
dx = 1.0;
xloc = 0:dx:(nx-1)*dx;
%Initial and boundary conditions
h0 = zeros(nx,1);
h0(1) = hL; h0(nx) = hR;
dh = (hL-hR)/((nx-1)*dx);
h0(2:nx-1) = hL-dh:-dh:hR+dh;

Q = zeros(1,nx);
Q(30) = -3e-2;
Q(60) = 3e-2;
Ss = 1e-5;
S = Ss*thick; %storativity
%Time
dt = 0.02;
tmax = 0.2;
ts = zeros(ceil(tmax/dt),1)+dt;
nt = length(ts);
% divide the domain into 3 parts
N=3;
% 
Ym = 0;
sigmaY = 1.0;
%integral scale
ind_scale = 10.0;
iscale = ind_scale;
%Set synthetic truth
truek = exp(randn(N,1)).*thick;
truekfield = zeros(1,nx);
truekfield(1:33) = truek(1);
truekfield(34:66) = truek(2);
truekfield(67:100) = truek(3);
%Solve truefield
figure(1); clf;
h = h0;
hold on;
for i=1:nt
    h = hsolver(h, dt, truekfield);
    plot(h, 'r', 'LineWidth', 1);
end
hold off;
xlabel('x (L)');
ylabel('Head (L)');

% 
% Use 2-order Hermite Polynomials
%
polyindx = [2 2 2];
P=3; 
[hmat xmat,h2mat] = genhmat(N,P, polyindx);
M = max(size(xmat));
kfield = zeros(1,nx);

%Assume k follows uniform distribution
%Define bounds for uniform distribution
lb = [0.1 0.1 0.01];
ub = [2   2   1];

figure(4), clf;
hold on;
%% do MC
Nsamples = 3e4;
itype = 1; %1=normal distribution, 2=uniform distribution
if (itype==2)
    umat = rand(Nsamples, 3);
elseif (itype==1)
    umat = randn(Nsamples, 3);
end
modout = zeros(Nsamples,nx);
for i=1:Nsamples
    if (itype==2) 
    kfield(1:33)  = lb(1)+(ub(1)-lb(1))*umat(i,1);
    kfield(34:66) = lb(2)+(ub(2)-lb(2))*umat(i,2);
    kfield(67:100)= lb(3)+(ub(3)-lb(3))*umat(i,3);
    elseif (itype==1)
    kfield(1:33)  = exp(Ym+sigmaY*umat(i,1));
    kfield(34:66) = exp(Ym+sigmaY*umat(i,2));
    kfield(67:100)= exp(Ym+sigmaY*umat(i,3));
    end
    kfield = kfield.*thick; 
    h = h0;

    for it=1:nt
        h = hsolver(h, dt, kfield);
    end
    modout(i,:) = h;
    plot(h);
end
meanHead = mean(modout,1);
varHead = var(modout, 0, 1);
sigmaHead = sqrt(varHead);
figure(7), clf;
subplot(2,1,1);
plot(meanHead);
subplot(2,1,2);
plot(varHead);
if (itype==2)
    save 'mc_uniform.mat' 'meanHead' 'varHead' 'sigmaHead'
elseif (itype==1)
    save 'mc_normal.mat' 'meanHead' 'varHead' 'sigmaHead'
end