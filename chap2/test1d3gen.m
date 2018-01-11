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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test1D
%Case 3
%Author: Alex Sun
%Date: 3/23/2005
%Project:Source ID
%Generates response matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [amat] = test1d3gen()
%Global variables
%t0 starting time
%t1 ending time
%c0 source concentration
%V velocity
%dp Dispersion coeff
%decay decay coefficient
%L total length of 1D domain
%dx block length
global V dp alphaL;
global L h dx;
global t0s t1s nPrd;
global obswells obstimes nwells ntimes nsamples;
%stores sample values corresponding to 1 unit input
samplesvals = zeros(1, nsamples);
amat=zeros(nsamples, nPrd);
sunits=zeros(1, nPrd);
for k=1:nPrd
    sunits(k)=1.0;
    samplevals=samplegen1d3(sunits);
    amat(:,k)=samplevals';
    sunits(k)=0.0;
end
