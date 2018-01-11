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
%For Example 6.4
%Purpose: Factor analysis using the dataset presented in Figure 6_11
%This plots Figure 6_11a&b
%Date: $20110513$
%Dependency:
%	Matlab statistics toolbox
%====================================================================
close all;
clear all;
corrmat = [
1.000	0.422	0.540	0.277	0.847	0.864	0.536	0.989	0.798
0.422	1.000	0.697	0.350	0.335	0.584	0.428	0.430	0.586
0.540	0.697	1.000	0.799	0.284	0.583	0.779	0.523	0.696
0.277	0.350	0.799	1.000	-.150	0.233	0.786	0.277	0.554
0.847	0.335	0.284	-.150	1.000	0.765	0.106	0.859	0.643
0.864	0.584	0.583	0.233	0.765	1.000	0.376	0.890	0.858
0.536	0.428	0.779	0.786	0.106	0.376	1.000	0.505	0.585
0.989	0.430	0.523	0.277	0.859	0.890	0.505	1.000	0.845
0.798	0.586	0.696	0.554	0.643	0.858	0.585	0.845	1.000
];
%% After dropping sulfate(i.e., 5th column), the matrix is full-rank
% and positive definite
colno=5;
if (colno==1)
    range = 2:9;
elseif (colno==9)
    range = 1:8;
else
    range = [1:colno-1 colno+1:9];
end
disp(colno);
cmat = corrmat(range,range);
%cmat must be symmetric & positive definite for the factoran to work
det(cmat)
%% Select 2 factors because there are only two eigenvalues>1
nfactors = 2;
[lambda,psi] = factoran(cmat,nfactors, ...
                            'xtype', 'covariance', ...
                            'rotate', 'varimax', ...
                            'nobs', 17);
% Figure 6.11(a)                        
figure(1), clf;
eigvals = sort(eig(cmat), 'descend');
plot(eigvals, '-o', 'LineWidth', 1.5);
xlabel('Eigenvalue Index');
ylabel('Eigenvalue');
figure(2), clf;
hold on;
box on;
plot(lambda(:,1), 'b', 'LineWidth', 1.5);
plot(lambda(:,2), 'r', 'LineWidth', 1.5);
ylabel('Factors')
set(gca,'XTickLabel',{'EC','TOC','BOD','COD','Sodium', 'T. Coli', 'TDS', ...
    'Total N'});
plot([1 5 7 8],lambda([1,5,7,8],1), 'bs','MarkerFaceColor','b');
plot([2 3 4 6],lambda([2 3 4 6],2), 'ro');
hold off;
%% Figure 6.11(b)

figure(3), clf;
h=bar(lambda, 'BaseValue', 0.5);
set(gca,'XTickLabel',{'EC','TOC','BOD','COD','Sodium', 'T. Coli', 'TDS', ...
    'Total N'});
set(get(h(1),'BaseLine'),'LineWidth',2,'LineStyle',':')
ylabel('Factor Loading');