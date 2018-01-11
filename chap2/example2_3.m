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
%Revision $20110817$
%for Example 2.3 and Example 2.4
%Dependency:
%   -test1d3gen, samplegen1d3
%Add uniform random noise
%Use breakpoint to view variables
%This example is used to illustrate
%(1) ill-conditioned & overdetermined system
%(2) The effect of measurement error
%====================================================================
function [] = example2_3
%====================================================================
global V alphaL dp;
global L h dx;
global t0s t1s nPrd;
global obswells obstimes nwells ntimes nsamples;
%Longitudinal dispersivity
alphaL=0.1;
%Advective velocity
V=1.0;
%Define release periods
t0s=[0  1 2 3 4 5 6];
t1s=[1  2 3 4 5 6 7];
c0s=[10 5 30 80 20 100 50];

%Number of release periods
nPrd = length(t0s);
%Define grid
%L-total length, dx=discretization
L=10.0;
h=1.0;
dx=0.2;
%Define observation wells
obswells=[1, 3, 5, 7, 9];
obstimes=[0.5,0.8,1.2,1.8,2.0,10.0];
nwells=max(size(obswells));
ntimes=max(size(obstimes));
nsamples=nwells*ntimes;
%Generate response matrix, A
trueA=test1d3gen;
[m,n]=size(trueA);
%Generate observation vector, b
trueb=samplegen1d3(c0s)';
truex=trueA\trueb
pinv(trueA)*trueb

%reset random generator seed to make results reproducible
rng = RandStream.create('mt19937ar','seed',69999);
RandStream.setGlobalStream(rng);

%case 1: 5% error
eamp=0.05;
b1= trueb.*(1-eamp+2*eamp*rand(size(trueb)));
x1 = trueA\b1
pinv(trueA)*b1
rmse1 = sqrt(sum((truex-x1).^2)./n)
%case 2: 10% error
eamp=0.1;
b2= trueb.*(1-eamp+2*eamp*rand(size(trueb)));
x2 = trueA\b2
pinv(trueA)*b2
rmse2 = sqrt(sum((truex-x2).^2)./n)
return
