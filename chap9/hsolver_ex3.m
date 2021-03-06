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
%Purpose: Defines the matrices required for extended KF example
%==========================================================================
%Parameter setup
function [y] = hsolver_ex3(y0,params)
% h0, result from previous step
% dt, time step
% T, transmissivity
global q S;
global nx dx;
global H0 qL;

Amat = zeros(nx, nx);
bvec = zeros(nx, 1);

h0 = y0(1:nx);
dt = params{1};
a1 = S/dt;
dx2 = 1/(dx*dx);
T = exp(y0(end));
a2 = T*dx2;

k = 1;
for i=1:nx
    if (i>1)
        Amat(k, k-1) = -a2;
    end
    if (i<nx)
        Amat(k, k) = a1+2*a2;
        Amat(k, k+1) = -a2;
    end
    bvec(k) = q(i);
    if (i==1)
        bvec(k) = bvec(k)+a2*H0;
    elseif (i==nx)
        Amat(k,k)=a1+a2;
        bvec(k) = bvec(k)-qL/dx;
    end
    k=k+1;
end
h = Amat\(a1*h0+bvec);
y = [h; y0(end)];