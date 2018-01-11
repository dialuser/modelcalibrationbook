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

%Figure 6.9
%Purpose: Illustrate relation between variogram and covariance
%Author: Alex Sun
%Date: $20130827$
%Dependency: requires mgstat library
%==========================================================================
addpath('./mGstat');
V1 = '1 Gau(1)';
figure(1), clf;
[sv, d] = semivar_synth(V1,0:0.05:1.5);
cv = 1-sv;
[AX, H1, H2] = plotyy(d, sv, d, cv);
set(H1,'LineWidth', 1.5);
set(H2,'LineWidth', 1.5);
set(H2,'LineStyle', '--');
set(get(AX(1),'Ylabel'),'String','\gamma(r)') 
set(get(AX(2),'Ylabel'),'String','C(r)') 
xlabel('Separation distance, r')
box on;