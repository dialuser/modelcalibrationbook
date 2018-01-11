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
%Figure 6.3
%Plot IDW illustration
%==========================================================================
function figure6_3
close all;
%reset random generator
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',12311));
%generate random data sites
data = randn(10, 1);
coord = rand(10, 2).*10;
%define fine grid for interpolation
x1 = 0:0.1:10;
y1 = x1;
[X1, Y1] = meshgrid(x1,y1);
[m,n] = size(X1);
%
%generate interpolation basis for different beta values
%
figure(1), clf;
subplot(1,3, 1);
hold on;
b = 1.0;
[Z1,W1,dist1] = idw(b);
%mesh(X1,Y1, Z1);
surf(X1,Y1,Z1, 'LineStyle', 'none');
plot3(coord(:,1), coord(:,2), data, 'ko', 'LineWidth', 1.5);
zlim([-2 1.5]);
title(['\beta=' num2str(b,'%2.1f')], 'fontsize', 12, 'fontname', 'arial');
viewsetting = [-32 18];
view(viewsetting);
grid on;
hold off;
subplot(1,3,2);
hold on;
b = 2;
[Z2,W2,dist2] = idw(b);
%mesh(X1,Y1, Z2);
surf(X1,Y1,Z1, 'LineStyle', 'none');
zlim([-2 1.5]);
plot3(coord(:,1), coord(:,2), data, 'ko', 'LineWidth', 1.5);
title(['\beta=' num2str(b,'%2.1f')], 'fontsize', 12, 'fontname', 'arial');
view(viewsetting);
grid on;
hold off;
subplot(1,3, 3);
hold on;
b = 4;
[Z3,W3,dist3] = idw(b);
%surfl(X1,Y1, Z3);
surf(X1,Y1,Z1, 'LineStyle', 'none');
shading interp;
plot3(coord(:,1), coord(:,2), data, 'ko', 'LineWidth', 1.5);
zlim([-2 1.5]);
view(viewsetting);
title(['\beta=' num2str(b,'%2.1f')], 'fontsize', 12, 'fontname', 'arial');
grid on;
hold off;
    function [Z,W11,dist11] = idw(b)
        %this function does actual idw
        Z = zeros(n,m);
        for i=1:m
            for j=1:n
                dist = sqrt((x1(i)-coord(:,1)).^2 + (y1(j)-coord(:,2)).^2);
                W = dist.^(-b);
                W = W./sum(W);
                Z(j,i) = sum(W.*data);
                if (i==40 && j==40)
                    W11=W;
                    dist11=dist;
                end
            end
        end

    end
end