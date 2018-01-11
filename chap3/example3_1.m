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
%Date: $20130616$
%For Example 3.1
%Dependency:
% 	Requires Matlab optimization toolbox
%Modified to use objective function norm^2
%Modified parameters and coupled equation
%====================================================================
function example3_1()
%Example 3_1: estimate parameters in Monod equations
%Define true parameters
%theta(1): mu_max, maximum substrate utilization rate per unit biomass
%theta(2): Ks, half saturation constant
%theta(3): Y, yield coefficient
%b: endogenous decay rate (assumed known)
thetat=zeros(3,1);
%the true values
thetat(1) = 0.15;
thetat(2) = 0.4;
thetat(3) = 0.2;
%define total simulation time
ttotal = [0 30];
initconc = [19 0.05];
b = 0.0042;
%make synthetic truth
theta= thetat;
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
% This is Figure 3.2
figure(1), clf;
subplot(1,2,1);
plot(T0,conc0(:,1),'-', Tobs, concObs(:,1), 'o', 'LineWidth', 1.5);
set(gca, 'FontName', 'Arial','FontSize', 12);
xlabel('Time (h)', 'fontsize',12, 'fontname', 'arial');
ylabel('Substrate Conc. (mg/L)', 'fontsize',12, 'fontname', 'arial');
legend('True', 'Observation', 'Location', 'Southwest');
legend boxoff;
subplot(1,2,2);
plot(T0,conc0(:,2),'-', Tobs, concObs(:,2), 'o', 'LineWidth', 1.5);
set(gca, 'FontName', 'Arial','FontSize', 12);
xlabel('Time (h)','fontsize',12, 'fontname', 'arial');
ylabel('Biomass Conc. (mg/L)', 'fontsize',12, 'fontname', 'arial');

% Prior information of parameters
theta0 = [0.1 0.5 0.1];
% weighting coefficient
alphaSet = [0.01:0.02:0.5];
objval1 = zeros(size(alphaSet));
objval2 = zeros(size(alphaSet));

for ik = 1:length(alphaSet)
    alpha = alphaSet(ik)
    x0 = theta0; %initial guess    
    thetaparam = fminunc(@objfun, x0);
    objval2(ik) = norm(thetaparam-theta0).^2;
    theta=thetaparam
    %this is used to calculate distance to the true solution
    norm(theta-thetat')
    [T,concSol] = ode45(@monod,ttotal,initconc);
    concModel = interp1(T, concSol, Tobs);
    objval1(ik) = norm(concModel-concObs).^2;
end

figure(2), clf;
%This is Figure 3.3 
plot(objval2,objval1, 'bx', 'LineWidth', 1.5);
xlabel('||\theta-\theta_{0}||^2', 'fontsize',12, 'fontname', 'arial');
ylabel('||u(\theta)-u_{obs}||^2', 'fontsize',12, 'fontname', 'arial');
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);

%Use the following to plot the final plot
%This corresponds to alpha=0.31
figure(3), clf;
theta=[  0.1453    0.1177    0.1885];
[T,concSol] = ode45(@monod,ttotal,initconc);
subplot(1,2,2);
plot(T0, conc0(:,2),'-', Tobs, concObs(:,2), 'o', T, concSol(:,2), '--', 'LineWidth', 1.5);
subplot(1,2,1);
plot(T0, conc0(:,1),'-', Tobs, concObs(:,1), 'o', T, concSol(:,1), '--', 'LineWidth', 1.5);
legend('True', 'Observation', 'Bi-criterion', 'Location', 'Southwest');
legend boxoff;
disp('successfully finished Example 3.1');
%% coupled ode
%% coupled ode
    function dc = monod(t,c)
    %c(1): substrate concentration
    %c(2): biomass concentration
    dc = zeros(2,1);
    mu = theta(1)*c(1)./(theta(2)+c(1));
    dc(1) = -mu/theta(3).*c(2);
    dc(2) = mu.*c(2)-b*c(2);
    end
%% objective func
    function [res] = objfun(x)
        theta=x;
        [T, concAll] = ode45(@monod,ttotal,initconc);
        %get observed data
        concModel = interp1(T, concAll, Tobs);
        res = norm(concObs-concModel).^2 + alpha*norm(x-theta0).^2;
    end
end