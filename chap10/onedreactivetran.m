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
%One-d transport model
function [res]=onedreactivetran(par)
%van Genutchen and Parker 1984
%R is retardation factor
%Flux boundary (v*c) at x=0, t>0
%dc/dx=0 at x=L
% domain length
global x t;
V =par(1); %darcy flux
theta = par(2); %porosity
V = V/theta;
alpha=par(3); %longitudinal dispersivity
De=par(4); %effective diffusion coeff
D = alpha*V+De;
R = par(5);

term1 = (R.*x-V.*t)./(2*sqrt(D*R*t));
term2 = (R.*x+V.*t)./(2*sqrt(D*R*t));
B =0.5*erfc(term1)+sqrt(V.^2.*t./(pi*D*R)).*exp(-term1.^2) ...
    -0.5*(1+V.*x/D+V.^2.*t./(D*R)).*exp(V.*x./D).*erfc(term2);

res = B;
