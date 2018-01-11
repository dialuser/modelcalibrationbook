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

%Author:Alex Sun
%Date: $20120118$
%For Example 3.3
%
%Dependency:
% 	Requires Matlab optimization toolbox
%====================================================================
function example3_3()
%Example 3_3: estimate parameters in Monod equations
%by solving a contrained optimization problem
%Define true parameters
%theta(1): mu_max, maximum substrate utilization rate per unit biomass
%theta(2): Ks, half saturation constant
%theta(3): Y, yield coefficient
%b: endogenous decay rate (assumed known)
thetat=zeros(3,1);
thetat(1) = 0.15;
thetat(2) = 0.4;
thetat(3) = 0.2;
theta= thetat;
ttotal = [0 30];
initconc = [19 0.05];
b = 0.0042;
%make synthetic truth
[T0,conc0] = ode45(@monod,ttotal,initconc);
%make "contaminated" observations
rng = RandStream.create('mt19937ar','seed',39999);
RandStream.setGlobalStream(rng);
Tobs = 2:2:30;
N = length(Tobs);
concObs = interp1(T0, conc0, Tobs);
concObs(:,1) = concObs(:,1) + 0.2*randn(N, 1);
concObs(:,2) = concObs(:,2) + 0.02*randn(N, 1);
% plot true and noisy concentrations
figure(1), clf;
subplot(1,2,1);
plot(T0,conc0(:,1),'-', Tobs, concObs(:,1), 'o', 'LineWidth', 1.5);
set(gca, 'FontName', 'Arial','FontSize', 12);
xlabel('Time (h)', 'fontsize',12, 'fontname', 'arial');
ylabel('Substrate Conc. (mg/L)', 'fontsize',12, 'fontname', 'arial');
legend('True', 'Observation');
legend boxoff;
subplot(1,2,2);
plot(T0,conc0(:,2),'-', Tobs, concObs(:,2), 'o', 'LineWidth', 1.5);
set(gca, 'FontName', 'Arial','FontSize', 12);
xlabel('Time (h)','fontsize',12, 'fontname', 'arial');
ylabel('Biomass Conc. (mg/L)', 'fontsize',12, 'fontname', 'arial');

%% Solve constrained optimization
% first set linear bounds for all three parameters
lb = [0.1 0.1 0.1];
ub = [1.0 1.0 0.5];
x0 = [0.1 0.5 0.1]; %initial guess
theta0 = [0.15 0.35 0.25]; 

alpha = 0.01;
options=optimset('Algorithm','active-set', 'TolFun', 1e-7);
[thetaSolution,fval,exitflag,output,lambda]=fmincon(@objfun,x0,[],[],[],[],lb,ub,@confun, options)
figure(3), clf;
theta=thetaSolution;
[T,concSol] = ode45(@monod,ttotal,initconc);
subplot(1,2,1);
plot(T0, conc0(:,2),'-', Tobs, concObs(:,2), 'o', T, concSol(:,2), '--', 'LineWidth', 1.5);
subplot(1,2,2);
plot(T0, conc0(:,1),'-', Tobs, concObs(:,1), 'o', T, concSol(:,1), '--', 'LineWidth', 1.5);

alpha = 0.2;
[thetaSolution2,fval,exitflag,output,lambda2]=fmincon(@objfun,x0,[],[],[],[],lb,ub,@confun, options)
figure(4), clf;
theta=thetaSolution2;
[T,concSol] = ode45(@monod,ttotal,initconc);
subplot(1,2,1);
plot(T0, conc0(:,2),'-', Tobs, concObs(:,2), 'o', T, concSol(:,2), '--', 'LineWidth', 1.5);
subplot(1,2,2);
plot(T0, conc0(:,1),'-', Tobs, concObs(:,1), 'o', T, concSol(:,1), '--', 'LineWidth', 1.5);

disp('successfully finished Example 3.3');



%% coupled ode
    function dc = monod(t,c)
    %c(1): substrate concentration
    %c(2): biomass concentration
    dc = zeros(2,1);
    mu = theta(1)*c(1)./(theta(2)+c(1));
    dc(1) = -mu/theta(3).*c(2);
    dc(2) = mu*c(2)-b*c(2);
    end
%% objective func
    function [res] = objfun(x)
        theta=x;
        [T, concAll] = ode45(@monod,ttotal,initconc);
        %get observed data
        concModel = interp1(T, concAll, Tobs);
        res = (norm(concObs-concModel)).^2;
    end
%% c fun
    function [c,ceq] = confun(x)
        c = norm(x-theta0)- alpha;
        ceq=[];
    end
end