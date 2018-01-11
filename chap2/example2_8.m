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

%Example 2.8
%Author: Alex Sun
%Purpose: For example 2.8: estimate van Genutchen parameters using
%different norms
%Dependency:
%Requires Matlab optimization toolbox
%====================================================================
%Generate obs data
function example2_8()
theta_s=0.403; %[]
theta_r=0.029; %[]
%Experimental data from 
%van Genuchten, M.Th. 1980. A closed-form equation for predicting the
%hydraulic conductivity of unsaturated soils. Soil Sci. Soc. Am. J.
%44:892–899.
data=[0.03 440;
      0.04 170;
      0.05 70;
      0.06 34;
      0.08 13;
      0.10 4.5;
      0.12 3.3;
      0.14 2.59;
      0.16 2.09;
      0.18 1.68;
      0.2  1.34;
      0.22 1.06;
      0.24 0.78;
      0.26 0.64;
      0.28 0.53;
      0.3  0.43;
      0.32 0.34;
      0.34 0.26;
      0.36 0.18;
      0.38 0.10;
      0.4  0.03;
      0.41 1e-2];
theta_data = data(:,1);
psi_data = data(:,2);
%Initial guess
x0=[2 1];
%Set options for fminunc
options=optimset('Display', 'iter', 'MaxFunEvals', 500);
solvertype = 2;
[xunc2]=fminunc(@myobj, x0 ,options)
solvertype = 1;
[xunc1]=fminunc(@myobj, x0 ,options)
solvertype = 'inf';
[xuncinf]=fminunc(@myobj, x0 ,options)
figure(1), clf;
theta_calibrated = vanGen(xunc1);
semilogy(theta_data, psi_data, 'o');
hold on;
semilogy(theta_calibrated, psi_data,'-', 'LineWidth', 1.5);
resval1= norm(theta_calibrated-theta_data)
theta_calibrated = vanGen(xunc2);
resval2= norm(theta_calibrated-theta_data)
semilogy(theta_calibrated, psi_data,'g-.', 'LineWidth', 1.5);
theta_calibrated = vanGen(xuncinf);
resvalinf= norm(theta_calibrated-theta_data)
semilogy(theta_calibrated, psi_data,'r--', 'LineWidth', 1.5);
hold off;
legend('Exp. Data', 'L_1', 'L_2', 'L_\infty');
legend boxoff;
ylabel('|\psi| [m]');
xlabel('Water content [-]');
% solvertype = 'KL';
% xunc4 =fminunc(@myobj, x0 ,options)
%% calculate objective function for fminunc
    function f=myobj(x)
        m = 1-1/x(2); 
        theta_e=(1+(x(1)*abs(psi_data)).^x(2)).^(-m);
        theta_model=theta_e.*(theta_s-theta_r)+theta_r;
        if (solvertype==2)
            f = norm(theta_data-theta_model);
        elseif (solvertype==1)
            f = norm(theta_data-theta_model, 1);
        elseif (strcmp(solvertype,'inf'))
            f = norm(theta_data-theta_model, inf);
        elseif (strcmp(solvertype,'KL'))
            f = sum(log(theta_model)-log(theta_data));
        end
    end
%% calculates matric potential for given VG parameter
    function f=vanGenP(x)
        s= (theta_data-theta_r)./(theta_s-theta_r);
        m = 1-1/x(2);
        f= 1/x(1).*(s.^(-1/m)-1).^(1/x(2));
    end
%% Generate effective saturation for given VG parameters and exp \psi data
    function f=vanGen(x)
        m = 1-1/x(2);
        f=(1+(x(1)*abs(psi_data)).^x(2)).^(-m);
        f= f.*(theta_s-theta_r)+theta_r;
    end
    end