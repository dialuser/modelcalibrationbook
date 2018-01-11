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
%Figure 6.5
%Shows Voronoi diagram and natural neighbor search
%Need to process in Illustrator to add the dotted lines inside the shaded
%area
%==========================================================================
%Generate random voronoi generators
%reset random generator
close all;
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',12320));
%generate random data sites
x = gallery('uniformdata',[1 20],0);
y = gallery('uniformdata',[1 20],1);

figure(1), clf;
subplot(1,2,1);
%Plot bounded voronoi diagram
voronoi(x,y); 
hold on;
plot(x, y, 'ko', 'markersize', 3);
hold off;
box on;
axis square;

subplot(1,2,2);
x=[x 0.55];
y=[y 0.55];
%Plot bounded voronoi diagram
voronoi(x,y); 

%This uses the method from Matlab voronoi manual page
for i = length(c) 
    if all(c{i}~=1)   % If at least one of the indices is 1, 
                      % then it is an open region and we can't 
                      % patch that.
       patch(v(c{i},1),v(c{i},2),i); % use color i.
    end
end
box on;
axis square
hold on;
plot(x(1:end-1), y(1:end-1), 'ko', 'markersize', 3);
plot(x(end), y(end), 'rx', 'markersize', 6, 'Linewidth', 1.5);
[v,c]=voronoin([x' y']); 
hold off;
