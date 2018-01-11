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
%Used by example6_5.m
%====================================================================
function fval = klchar(w)
%The characteristic function used for determining eigenvalues of 
%KL expansion
global iscale L;
eta = iscale/L;
w1 = w*L;
fval = (eta^2*w1.^2-1).*sin(w1)-2*eta.*w1.*cos(w1);
end

