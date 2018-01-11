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
%Revision $20130718$
%Figure 6.6
%Shows the Gaussian and Quadric windows
%Dependency:
%	Requires Matlab statistics toolbox
%====================================================================
figure(1), clf;
subplot(1,2,1);
%Gaussian
x=-3:0.1:3;
y = normpdf(x, 0, 1);
plot(x, y, 'LineWidth', 1.5)
axis square
subplot(1,2,2);
%quadric
beta = 2;
y = sqrt(1-(x./beta).^2);
ind=find(real(y)>0);
plot(x(ind), y(ind), 'LineWidth', 1.5)
axis square