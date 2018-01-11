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
%
%Author: Alex Sun
%Example 10.2
%Dependencies: lhsdesign, xvecconvert (included), onedreactivetran(included)
%Purpose: Demonstrate GSA algorithm using a simple 1D transport problem
%Color post-processed in Illustrator
%%
clear all; close all;
%make the results repeatable
rng('default');
rng(1);
global x t;
%param is parameter vector consisting of
%[darcy flux, porosity, longitudinal dispersivity, effective diffusion, retardation]
%All parameters are triangular distributionsf
parm  =[1e-5,   0.2,  5,   2e-9, 7];
parmlb=[1e-6,  0.10,  1, 1.8e-9, 5];
parmub=[5e-5,  0.25, 10, 2.2e-9, 15];
day2s=86400;
x=[0.01:0.01:1 2:20]; %[m]
t=10*day2s; %[s]
%This plots Figure 10.4 of the book
figure(1), clf;
hold on;
%show monte carlo effect using LHS
nMC = 1000;
MCPoints = lhsdesign(nMC, length(parm));
AllRes = zeros(nMC, length(x));
for ii = 1:nMC
   % convert random numbers to parameter values
   xvec = xvecconvert(parmlb, parmub, parm, MCPoints(ii,:));
   AllRes(ii,:) = onedreactivetran(xvec);
end
plot(x, AllRes, 'color', [0.7,0.7,0.7]);
conc = onedreactivetran(parm);
plot(x, conc, 'LineWidth', 2.0);
xlabel('x (m)');
ylabel('C/C_0');
ylim([0 1]);
hold off;
box on;

% Define different distances
allVi = [];
allVTi = [];
xval=[0.1 5 10 20];

for x = xval
    disp(['x= ' num2str(x)]);
    sumyA = 0;
    sumyA2 = 0;
    sumyAB = 0;

    k=length(parm);
    Vi = zeros(k, 1); 
    VTi = zeros(k, 1);

    %number of model runs
    N=2^13;
    TAll = lhsdesign(N, 2*k);
    for i=1:N
       %generation of uniform samples in (0;1)
       %using pseduorandom sequences <note: this requires LPTAU51 by
       %Satelli
       %T = LPTAU51(i,2*k);
       T = TAll(i,:);
       %PREPARATION OF THE RADIAL SAMPLE MATRIX X(k+2,k) 
       A=T(1:k);
       B=T(k+1:2*k);
       % convert random numbers to parameter values
       xvec = xvecconvert(parmlb, parmub, parm, A);
       yA = onedreactivetran(xvec);
       sumyA  = sumyA + yA;
       sumyA2 = sumyA2+ yA^2;
       xvec = xvecconvert(parmlb, parmub, parm, B);
       yB = onedreactivetran(xvec);

       sumyAB  = sumyAB + yA*yB;
       for jj=1:k
            Ab=B;
            Ab(jj)=A(jj);
            xvec = xvecconvert(parmlb, parmub, parm, Ab);
            yAb=onedreactivetran(xvec);
            Vi(jj) = Vi(jj) +yA.*yAb;
            VTi(jj)= VTi(jj)+yB.*yAb;
       end
    end %N
    N1=N-1;
    sumyA = sumyA./N;
    sumyAB = sumyAB./N1;

    f2 = sumyA.^2;
    sumyA2 = sumyA2./N1; 
    Vy=sumyA2 -f2;
    for jj=1:k
       Vi(jj)=(Vi(jj)./N1-sumyAB)./Vy;
       VTi(jj)=1-(VTi(jj)./N1-sumyAB)./Vy;
    end
    [Vi VTi]
    allVi = [allVi Vi];
    allVTi = [allVTi VTi];    
end
%
figure(2), clf;
subplot(1,3,1);
barh([allVi(:,1) allVTi(:,1)]);
xlim([0 1]);
set(gca, 'YTickLabel', {'v', 'phi', 'lambda', 'D_e', 'R'});
set(gca, 'FontName', 'Arial');
legend('First order', 'Total effect');
set(gca, 'FontName', 'Arial');
title('x=0.1 m');

subplot(1,3,2);
barh([allVi(:,3) allVTi(:,3)]);
xlim([0 1]);
set(gca, 'YTickLabel', {'v', 'phi', 'lambda', 'D_e', 'R'});
title('x=10 m');

subplot(1,3,3);
barh([allVi(:,4) allVTi(:,4)]);
xlim([0 1]);
set(gca, 'YTickLabel', {'v', 'phi', 'lambda', 'D_e', 'R'});
title('x=20 m');
