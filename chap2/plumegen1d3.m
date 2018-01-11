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
function [res]= plumegen1d3(t, prdconcs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test1D
%Case 3
%Author: Alex Sun
%Date: 3/23/2005
%Project:Source ID
%Generates a synthetic plume
%Analytical solution taken from 
%C.W. Fetter, Contaminant hydrogeology, 1992
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Params
%t, current time
%prdconcs, source concentration during each release period
%
%Global variables
%t0 starting time
%t1 ending time
%c0 source concentration
%V velocity
%dp Dispersion coeff
%L total length of 1D domain
%dx block length
%alphaL, longitudinal dispersion coefficient
global V alphaL dp;
global L h dx;
global t0s t1s nPrd;
dp = alphaL*V;
nx = ceil(L/dx);
ddx=[1:nx]*dx;
conc = zeros(1, nx);
for j=1:nPrd
    %use vector operation here
    coeff2=ddx.*V./dp;
    if (t>t0s(j))
          coeff=2.0*sqrt(dp*(t-t0s(j)));
          coeff3=V*(t-t0s(j));
          conc=conc+0.5*prdconcs(j).*( ...
          erfc((ddx-coeff3)./coeff)  ...
          +exp(coeff2).*erfc((ddx+coeff3)./coeff));
    end
    %Use mirror principle, to add a source of negative strength
    %to offset the positive source in the above
    if (t>t1s(j))
          coeff=2.0*sqrt(dp*(t-t1s(j)));
          coeff3=V*(t-t1s(j));
          conc=conc-0.5*prdconcs(j).*( ...
          erfc((ddx-coeff3)./coeff)  ...
          +exp(coeff2).*erfc((ddx+coeff3)./coeff));       
    end
end
%form return param
res=[ddx', conc'];