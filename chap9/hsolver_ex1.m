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
%==========================================================================
function [h Amat bvec cvec] = hsolver_ex1(h0, dt, T)
%Author: Alex Sun
%Purpose: Defines the matrices required for KF example
%         The left-hand-side boundary condition is asssumed unknown in this
%         example
%INPUT:
%   h0, result from previous step (i.e., k-th step)
%   dt, time step
%   T, transmissivity
%OUTPUT:
%   h, result at (k+1)-th step
%   Amat, the coefficient matrix
%   bvec, vector containing forcing/control variables (known)
%   cvec, vector containing forcing/control variables (unknown)
global q S;
global nx dx;
global H0 HL;
global a1;

Amat = zeros(nx, nx);
bvec = zeros(nx, 1);
cvec = eye(nx, 1);
a1 = S/dt;
a2 = T/(dx*dx);

k = 1;
for i=1:nx
    Amat(k, k) = a1+2*a2;
    if (i>1)
        Amat(k, k-1) = -a2;
    end
    if (i<nx)
        Amat(k, k+1) = -a2;
    end
    bvec(k) = q(i);
    if (i==nx)
        bvec(k) = bvec(k)+a2*HL;
    end
    k=k+1;
end

cvec(1) = a2;
h = Amat\(a1*h0+bvec+cvec*H0);