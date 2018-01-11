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
%Transform random variables back to their actual values
%from uniform random variates
%All variates are triangular
function [ xvec ] = xvecconvert(xlb, xub, xmode, rvec)
    xvec = zeros(size(xlb));
    Fc = (xmode-xlb)./(xub-xlb);
    for i=1:length(xvec)
        if (rvec(i)<Fc(i))
            xvec(i) = xlb(i)+sqrt(rvec(i)*(xub(i)-xlb(i))*(xmode(i)-xlb(i)));
        else
            xvec(i) = xub(i)-sqrt((1-rvec(i))*(xub(i)-xlb(i))*(xub(i)-xmode(i)));
        end
    end
end

