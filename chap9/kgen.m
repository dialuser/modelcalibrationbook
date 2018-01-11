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
%Date: $20080308$
%Generate correlated random field
%==========================================================================
function [karr] = kgen(Ym, sigmaY, nx, xloc, iscale, nrz, nongau, ...
    imethod)
if (~exist('nongau', 'var'))
    covmat = gencov(nx, xloc, iscale, sigmaY);
    [V, D] = eig(covmat);
    % Generate random realizations for Y
    karr = Ym + V*sqrt(D)*randn(nx, nrz);
    %convert to K
    karr = exp(karr);
else
    %the first component of Ym is indicator
    %the second is Y1
    %the third is Y2    
    if (length(Ym) ~= 3) 
        error('Invalid input variable dimension');
    end

    if (imethod == 1)
        %generate continuous standard normal variates
        %using exponential model
        covmat = gencov(nx, xloc, iscale(1), 1.0);
        [V, D] = eig(covmat);
        % Generate random realizations
        karr1 = V*sqrt(D)*randn(nx, nrz);
        % Do truncation at invnormal(P = 0.5)
        cutoff = icdf('Normal', Ym(1), 0, 1);    
        indx2 = find(karr1>cutoff);
        % Simulate kfield 1
        covmat = gencov(nx, xloc, iscale(2), sigmaY(2));
        [V, D] = eig(covmat);
        karr1 = Ym(2) + V*sqrt(D)*randn(nx, nrz);    
        covmat = gencov(nx, xloc, iscale(3), sigmaY(3));
        [V, D] = eig(covmat);
        karr2 = Ym(3) + V*sqrt(D)*randn(nx, nrz);    
        karr1(indx2) = karr2(indx2);
        karr = exp(karr1);
    else
        % Simulate kfield 1
        covmat = gencov(nx, xloc, iscale(2), sigmaY(2));
        [V, D] = eig(covmat);
        karr1 = Ym(2) + V*sqrt(D)*randn(nx, nrz);    
        covmat = gencov(nx, xloc, iscale(3), sigmaY(3));
        [V, D] = eig(covmat);
        karr2 = Ym(3) + V*sqrt(D)*randn(nx, nrz);    
        karr = (exp(karr1)+exp(karr2));
    end
end
return;
%==========================================================================
function [covmat] = gencov(nx, xloc, iscale, sigmaY)

covmat = zeros(nx,nx);
for i=1:nx
    for j=i:nx
        covmat(i,j) = abs(xloc(j)-xloc(i));
        covmat(j,i) = covmat(i,j);
    end
end
covmat = -covmat./iscale;
covmat = sigmaY^2.*exp(covmat);


