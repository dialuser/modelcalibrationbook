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

%Figure 6.10
%Purpose: Plots three variogram models
%Author: Alex Sun
%Date: $20110420$
%Note: spherical variogram integral is set to 3 so that it matches 
%      the exponential and gaussian
%Dependency: mgstat (included)
%Rev:
%Date: $20130828$
%   remove Matern from the plot
%
%==========================================================================
addpath('./mGstat');
ismatern=0;
V1 = '1 Sph(1)';
V2 = '1 Gau(1)';
V3 = '1 Exp(1)';
if (ismatern)
    V4.type = 'Matern';
    V4.par1 = 1;
    V4.par2 = 1;
    V4.par3 = 1.5;
end
ub = 2.0;
figure(1), clf;
hold on;
[sv, d] = semivar_synth(V1,[0:0.1:ub]);
plot(d, sv, 'b', 'LineWidth', 1.7);
[sv, d] = semivar_synth(V2,[0:0.1:ub]);
plot(d, sv, 'r--', 'LineWidth', 1.7);
[sv, d] = semivar_synth(V3,[0:0.1:ub]);
plot(d, sv, 'g:', 'LineWidth', 1.7);
if (ismatern)
    [sv, d] = semivar_synth(V4,[0:0.1:ub]);
    plot(d, sv, 'k-.', 'LineWidth', 1.7);
end
hold off;
if (ismatern)
    legend('Spherical', 'Gaussian', 'Exponential', 'Matern', 'Location', 'SouthEast');
else
    legend('Spherical', 'Gaussian', 'Exponential', 'Location', 'best');
end
legend boxoff;
xlabel('r');
ylabel('\gamma(r)')
box on;