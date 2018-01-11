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
%Purpose: Solve the characteristic equation to get roots
%For Example 6.5
%Last update:
%$01122014$
%====================================================================
%% Part(a): plots the first 4 eigenfunctions
close all; clear all;
global iscale L;
iscale = 1;
L = 10;
sigma2=1.0;
x = 0:0.1:L;
frange = klchar(x);
figure(1), clf;
plot(x, frange);
nroots = 50;
% Get non-trivial roots only
intervs = [0.3 0.5 0.8 1.0 1.3 1.6 1.9 2.2 2.5 2.8 3.2 3.5 3.8 4.1 4.4 4.7 5.0 5.3 5.6 6.0]; 
%
w = findallzeros(@(x)klchar(x),0.1, nroots);
w = w';

lambda=2*iscale*sigma2./(iscale^2*w.^2+1);

x=[0:0.1:L];
if (nroots==4)
    figure(2), clf;
    box on;
    hold on;
    styles={'-', '--', '.-', ':'};
end
fmat = zeros(nroots, length(x));
for i=1:nroots
    fi = (iscale.*w(i)*cos(w(i)*x)+sin(w(i)*x))/...
        sqrt((iscale^2*w(i)^2+1)*L/2+iscale);
    if (nroots==4)
        plot(x, fi, styles{i}, 'LineWidth', 1.2);
    end
    fmat(i,:)=fi;
end
if (nroots==4)
    legend('1st', '2nd', '3rd', '4th', 'Location', 'SouthWest');
    legend boxoff;
    hold off;
    xlabel('Length, x');
    ylabel('Eigenfunction, \psi_i(x)');
end
%% Part B: covariance 
covfunc = sigma2.*exp(-x./iscale);
% 
% Get the approximation of the covariance
fvec = fmat(:,1);
pmat = repmat(fvec, 1, length(x)).*fmat.*repmat(lambda', 1, length(x));
p = sum(pmat, 1);
% Save the results for plotting
filename = sprintf('kleroot%d', nroots);
save(filename, 'p');

if (exist('kleroot4.mat', 'file') && ...
    exist('kleroot20.mat', 'file') && ...
    exist('kleroot50.mat', 'file'))

    load('kleroot4.mat');
    res4 = p;
    load('kleroot10.mat');
    res10 = p;
    load('kleroot50.mat');
    res50 = p;
    figure(3), clf;
    %
    hold on;
    plot(x, covfunc, 'LineWidth', 2.0);
    plot(x, res4, 'r--', 'LineWidth', 1.5);
    plot(x, res10, 'r-.', 'LineWidth', 1.5);
    plot(x, res50, 'r:', 'LineWidth', 1.5);
    ylim([-0.5 1]);
    xlabel('Separation distance, r=|x_1-x_2|');
    ylabel('Cov(r)');
    box on;
    hold off;
    legend('True', '4-term', ...
           '10-term', '50-term', 'Location', 'Northeast')
    legend boxoff;
    % Approximation Error
    rmseerr1 = sqrt(mean((covfunc-res4).^2))
    rmseerr2 = sqrt(mean((covfunc-res10).^2))
    rmseerr3 = sqrt(mean((covfunc-res50).^2))
else
    error('Need to obtain results for nroots=4, 10, and 50, respectively, before the results can be plotted');
end
return
%% Part b: plots the approximation (not used)
RandStream.setDefaultStream ... 
     (RandStream('mt19937ar','seed',123374));
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;

%% Generate a high-order random function using length(x) terms
% 
% Generate covariance matrix
%
nnodes = length(x);
covmat = zeros(nnodes, nnodes);
xloc = x;
for i=1:nnodes
    for j=i:nnodes
        p1 = xloc(i);
        p2 = xloc(j);
        covmat(i,j) = norm(p1-p2);
        covmat(j,i) = covmat(i,j);
    end
end
covmat = -covmat./iscale;
covmat = exp(covmat);
%
% Eigenvalue decomp
%
[V, D] = eig(covmat);
%
% Generate random realizations
%
zvec = sqrt(D)*randn(nnodes,1);
xvec = V*zvec;
%
% KL expansion using the first 4 random variables
%
rvec = randn(4,1); 
z1 = repmat(rvec, 1, length(x));
pmat = fmat.*repmat(sqrt(lambda'),1, length(x)).*z1;
p = sum(pmat,1);
%
% Plot comparison
%
figure(4), clf;
hold on;
plot(x, xvec);
xlabel('Length, x');
ylabel('Random function,\Phi(x,\omega)');
plot(x, p, '--', 'LineWidth', 1.5); 
hold off;
legend('Original', 'KL Expansion', 'Location', 'Southeast');
legend boxoff;
