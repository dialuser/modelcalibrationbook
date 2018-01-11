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
%Figure 6.4
%Shows Voronoi diagram
%To generate a bounded voronoi diagram,
% define a triangle with 3 points such that the resulting triangle 
% encompass all other points. Make sure these 3 points are
% pretty far away from all the other points. Then add these 3 points to
% the Voronoi algorithm and use axis to focus on the original set of
% points.
%==========================================================================
close all;
%Generate random voronoi generators
%reset random generator
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',16313));
%generate random data sites
coord = rand(100, 2).*10;
x = coord(:,1), y=coord(:,2);
%as08222013, add three points to encompass all points
x=[x; -15; 25; 5]; y = [y; 0; 0; 15];
coord = [x y];

figure(1), clf;
subplot(1,2,1);
%Plot bounded voronoi diagram
[v,c]=voronoin(coord); 
%This uses the method from Matlab voronoi manual page
for i = 1:length(c) 
    if all(c{i}~=1)   % If at least one of the indices is 1, 
                      % then it is an open region and we can't 
                      % patch that.
        patch(v(c{i},1),v(c{i},2),i); % use color i.
    end
end
xlim([0 10]); ylim([0 10]);
xlabel('x'); ylabel('y');
set(gca, 'XTick', 0:2:10);
set(gca, 'YTick', 0:2:10);
hold on;
plot(x, y, 'wo', 'markersize', 3, 'markeredgecolor', 'k');
hold off;
box on;

subplot(1,2,2);
%plot Delaunay dual diagram
tri = delaunay(x,y);
triplot(tri,x,y);
hold on;
plot(x, y, 'wo', 'markersize', 3, 'markeredgecolor', 'k');
hold off;
xlim([0 10]); ylim([0 10]);
xlabel('x'); ylabel('y');
set(gca, 'XTick', 0:2:10);
set(gca, 'YTick', 0:2:10);
