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
%For Example 2.5 in Chapter 2
%Dependency:
%   -test1d3gen, samplegen1d3
%Add uniform random noise
%This example is used to illustrate Truncated SVD solution
%====================================================================
function [Lmat Vmat residual trueA] = example2_5
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
%Perform 
[U S V] = svd(trueA);
sigvals = diag(S);
plot(sigvals, 'o-', 'LineWidth', 1.5);
xlabel('Index');
ylabel('Singular Values');
condnum = sigvals(1)/sigvals(7);
truncatedcondnum = sigvals(1)/sigvals(6);
truncA = V(:,1:6)*diag(1./sigvals(1:6))*U(:,1:6)';
truncatedsvd=truncA*trueb
rmse = sqrt(sum((truex-truncatedsvd).^2)./n)

%add obs error
%reset random generator seed to make results reproducible
rng = RandStream.create('mt19937ar','seed',69999);
RandStream.setGlobalStream(rng);

%case 1: 5% error
eamp=0.05;
b1= trueb.*(1-eamp+2*eamp*rand(size(trueb)));
x1 = truncA*b1
rmse1 = sqrt(sum((truex-x1).^2)./n)

%case 2: 10% error
eamp=0.1;
b2= trueb.*(1-eamp+2*eamp*rand(size(trueb)));
x2 = truncA*b2
rmse2 = sqrt(sum((truex-x2).^2)./n)

return
