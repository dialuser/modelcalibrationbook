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
%Solves 1-d transient groundwater flow equation for one time step
%
%==========================================================================
%Parameter setup
function [h] = hsolver(h0, dt, T);
% h0, initial condition
% dt, time step
% T, transmissivity
global Q S;
global nx dx;

%ibound = 1, const head both ends
%ibound = 2, const flux on right end

ibound = 1;
if (ibound == 1)
    Amat = zeros(nx-2, nx-2);
    bvec = zeros(nx-2, 1);
    h = h0;

    b = dt/S;
    c = b/dx^2;
    k = 1;
    for i=2:nx-1    
        Amat(k, k) = 1+c*(T(i)+T(i-1));
        if (i>2)
            Amat(k, k-1) = -c*T(i-1);
        end
        if (i<nx-1)
            Amat(k, k+1) = -c*T(i);
        end
        bvec(k) = h(i)+Q(i)*b;
        if (i==2)
            bvec(k) = bvec(k)+c*T(i-1)*h(i-1);
        elseif (i==nx-1)
            bvec(k) = bvec(k)+c*T(i)*h(i+1);
        end
        k=k+1;
    end
    h(2:nx-1) = Amat\bvec;    
elseif (ibound == 2)
    
end


