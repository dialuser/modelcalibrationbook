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

%Example 6.3, Figure 6.11
%Variogram Example 
%Data was downloaded from 
% http://www.ittvis.com/language/en-us/productsservices/envi/tutorials.aspx
% METADATA
% 640 cols, 400 rows, 6 bands
%Author: Alex Sun
%Date: $20110425$
%REV: $20130828$
%Dependency: requires mgstat
%==========================================================================
%Load ENVI multiband image
close all;
addpath('./mGstat');

nx=640;
ny=400;
x=multibandread('can_tmr.img', ...
    [ny,nx,6], 'int8', 0, 'bsq', 'ieee-le');
img = x(150:350, 150:350, 1);
figure(1), clf;
imagesc(img);
axis equal;
box on;
axis off;
% make offset and trimmed copies of image
figure(2), clf;
% define range
rng = 1:2:80;
gamma = zeros(size(rng));
i=1;
v = length(rng);

xoffset = 0;
if xoffset < 0  % difference is symmetric so can force xoffset positive
    xoffset = -xoffset;
    yoffset = -yoffset;
end

for yoffset = rng
    if yoffset > 0
        imga = img(1+yoffset:end, 1+xoffset:end);
        imgb = img(1:end-yoffset, 1:end-xoffset);
    else
        imga = img(1:end+yoffset, 1+xoffset:end);
        imgb = img(1-yoffset:end, 1:end-xoffset);
    end

    d = imga - imgb;
    v(i) = 0.5*mean(d(:).^2);
    i=i+1;
end
figure(2), clf;
hold on;
plot(rng, v, 'ro', 'LineWidth', 1.5);
xlabel('h');
ylabel('\gamma(h)');

yoffset = 0;
v1 = v;
i=1;
for xoffset = rng
    if yoffset > 0
        imga = img(1+yoffset:end, 1+xoffset:end);
        imgb = img(1:end-yoffset, 1:end-xoffset);
    else
        imga = img(1:end+yoffset, 1+xoffset:end);
        imgb = img(1-yoffset:end, 1:end-xoffset);
    end

    d = imga - imgb;
    v1(i) = 0.5*mean(d(:).^2);
    i=i+1;
end
[sv, d]=semivar_synth('3 Nug(0)+ 30 Exp (15)', rng);
plot(d, sv, 'r--', 'LineWidth', 1.5);

plot(rng, v1, 'b+','LineWidth', 1.5);
[sv, d]=semivar_synth('3 Nug(0)+ 34 Exp (12)', rng);
plot(d, sv, 'g-.', 'LineWidth', 1.5);

box on;
legend('N-S raw', 'N-S fitted', 'E-W raw', 'E-W fitted', 'Location', 'SouthEast');
legend boxoff;
hold off;