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
%Example 10.7 m /least squares approach
%Purpose: One-dimensional PCM problem
%Dependency: 
% hsolver.m, genhvec.m genhmat.m, 
% example10_7_PCM_MC.m
%Change log:
%   07/16/2012
%   change for least squares, sample more than M points
%   Figure 10.7
%==========================================================================
function []=example10_7_LS()
%Global Parameters
%S storativity
%Ss specific storage
%Q injection/pumping rate
%hL domain length
%dx cell length
close all;
global Q S;
global nx dx;
global ts;
global h0;
%
addpath('../chap8/sampling');
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',123375));

nx = 100;
hL = 100;
thick = 10;
hR = 90;
dx = 1.0;
xloc = 0:dx:(nx-1)*dx;
%Initial and boundary conditions
h0 = zeros(nx,1);
h0(1) = hL; h0(nx) = hR;
dh = (hL-hR)/((nx-1)*dx);
h0(2:nx-1) = hL-dh:-dh:hR+dh;

Q = zeros(1,nx);
Q(30) = -3e-2;
Q(60) = 3e-2;
Ss = 1e-5;
S = Ss*thick; 
%Time
dt = 0.02;
tmax = 0.2;
ts = zeros(ceil(tmax/dt),1)+dt;
nt = length(ts);
% Divide the domain into 3 zones
% The hydraulic conductivity in each zone is a random variable
% N, random dimension
N=3;
%Assume k follows uniform distribution
%Define bounds for uniform distribution
lb = [0.1 0.1 0.01];
ub = [2 2 1];
% 
% Choice of polynomials for each random dimension
% (1=Hermite, 2=uniform)
%
%% Case 1 Legendre
polyindx = [2 2 2];

P=2;
[~, varH1, numRZ1] = generateResults(P);

P=3;
[~, varH2, numRZ2] = generateResults(P);

P=4;
[~, varH3, numRZ3] = generateResults(P);
%Load MC results, generated separately using example10_10_PCM_MC
mcRes = load('mc_uniform.mat');

%
% Plotting
%
figure(1), clf;
hold on;
plot(varH1, 'k-',   'LineWidth', 2.0);
plot(varH2, 'b--',   'LineWidth', 2.0);
plot(varH3, 'r:',   'LineWidth', 2.0);
plot(mcRes.varHead, 'o', 'MarkerSize', 6);
xlabel('x (L)');
ylabel('\sigma^2_h (L^2)');
box on;
ylim([0 5]);
hold off;
str1 = ['PCM P=2, N_p=' num2str(numRZ1)];
str2 = ['    P=3, N_p=' num2str(numRZ2)];
str3 = ['    P=4, N_p=' num2str(numRZ3)];

legend(str1, str2, str3, 'Monte Carlo', 'Location', 'NorthWest');
%
axis tight

%%

function [meanH, varH, numberRuns] = generateResults(porder)
    %input: 
    %porder: order of PC expansion
    %M: number of terms in PC expansion
    M = factorial(N+porder)/(factorial(N)*factorial(porder));
    if (porder==2)
        nrz = 6*M;
    elseif (porder==3)
        nrz = 4*M;
    else
        nrz = 80;
    end
    %allSamples = lhsu(zeros(1,N)-1, zeros(1,N)+1, nrz);    
    allSamples = -1+2*rand(nrz,3);
    hmat = zeros(nrz, M);
    h2mat = zeros(nrz, M);
    
    modout = zeros(nrz,nx);
    
    for ir = 1:nrz
        [hvec h2vec] = genhvec(N, M, porder, allSamples(ir,:), polyindx);
        hmat(ir,:) = hvec;
        h2mat(ir,:) = h2vec;
        kfield = zeros(1,nx);

        if (polyindx(1)==1)
            kfield(1:33) = lb(1)+(ub(1)-lb(1))* ...
                            (0.5+0.5*erf(sqrt(2)/2*allSamples(ir,1)));
        elseif (polyindx(1)==2)
            kfield(1:33) = lb(1)+0.5*(ub(1)-lb(1))*(allSamples(ir,1)+1);
        end
        if (polyindx(2)==1)
            kfield(34:66) = lb(2)+(ub(2)-lb(2))* ...
                            (0.5+0.5*erf(sqrt(2)/2*allSamples(ir,2)));
        elseif (polyindx(2)==2)
            kfield(34:66) = lb(2)+0.5*(ub(2)-lb(2))*(allSamples(ir,2)+1);
        end
        if (polyindx(3)==1)
            kfield(67:100) = lb(3)+(ub(3)-lb(3))* ...
                            (0.5+0.5*erf(sqrt(2)/2*allSamples(ir,3)));
        elseif(polyindx(3)==2)
            kfield(67:100) = lb(3)+0.5*(ub(3)-lb(3))*(allSamples(ir,3)+1);
        end

        kfield = kfield.*thick;                
        h = h0;

        for it=1:nt
            h = hsolver(h, dt, kfield);
        end
        modout(ir,:) = h;
    end
    coeff = hmat\modout;
    save('coeff2.mat', 'coeff');
    %mean
    meanH = coeff(1,:);
    %Variance
    varH = zeros(1,nx);
    h2vec = h2mat(1,:);
    for j=2:M
        varH = varH + coeff(j,:).^2.*h2vec(j);
    end
    numberRuns = nrz; 
end
end

