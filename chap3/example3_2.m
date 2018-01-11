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

% Example 3.3
% purpose: use phillips function to illustrate Tikhonov regularization
%
%Dependency:
% 	Requires regularization toolbox (included)
%====================================================================
addpath('./regularization_matlab');
rng = RandStream.create('mt19937ar','seed',69999);
RandStream.setGlobalStream(rng);
[A,b,x] = phillips(64);
[U,s,V] = csvd(A);
b0 = b+1e-3*randn(size(b));

figure(1), clf;
subplot(1,2,1);
set(gca, 'FontName', 'Arial','FontSize', 12);
semilogy(s, '-d', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
hold on;
semilogy(abs(U'*b), 'rx', 'LineWidth', 1.5);
semilogy(abs(U'*b0), 'go', 'LineWidth', 1.5);
hold off;
%Use symbol d to label the curves to be consistent with the text
legend('s_i', '|U^T_i d/s_i|','|U^T_i d0/s_i|', 'Location', 'NorthEast')
ylim([1e-5 1e3]);
xlim([1 64]);
subplot(1,2,2);
set(gca, 'FontName', 'Arial','FontSize', 12);
[reg_corner,rho,eta,reg_param] = l_curve(U,s,b0);
rng = 1:5:length(rho);
plot_lc(rho(rng), eta(rng), '-o');

[x_Tik rhoval etaval]= tikhonov(U,s,V,b0, 2e-3)
figure(2), clf;
plot(x, 'r-');
hold on;
plot(x, 'g--');
hold off;