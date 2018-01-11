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

%Example 4.4
%Author Alex Sun
%Illustrates the results of Example 4.4 (Figure 4.2)
%Dependency:
%Requires Matlab statistics toolbox
%====================================================================

function example4_4()
figure(1), clf;
handle = subplot(1,2,1);
subplot(1,2,1);
plotpdfs(handle, [20,10,18.7,3.5,18.84,3.30],1)
box on;
legend('p_0', 'Likelihood', 'p_{*}');
legend boxoff;
xlabel('Conc (mg/L)');
ylabel('PDF');
title('(a)');
handle = subplot(1,2,2);
plotpdfs(handle, [20,10, 30,15, 23.08, 8.32],2)
box on;
xlabel('Conc (mg/L)');
ylabel('PDF');
title('(b)');

function plotpdfs(handle, param, num)
hold on;
mu_0=param(1);
sigma_0=param(2);
mu_1=param(3);
sigma_1=param(4);
mu_2=param(5);
sigma_2=param(6);
x=-30:80;
priorpdf = pdf('norm', x, mu_0, sigma_0)
plot(handle, x, priorpdf, 'b', 'LineWidth', 1.5);
likefun = pdf('norm', x, mu_1, sigma_1)
plot(handle,x, likefun, 'b-.', 'LineWidth', 1.5);
posteriorpdf = pdf('norm', x, mu_2, sigma_2)
plot(handle,x, posteriorpdf, 'b:', 'LineWidth', 1.5);
hold off;