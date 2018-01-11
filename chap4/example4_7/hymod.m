function [objval output] = hymod(parVec)
%Hymod.m
%Implementation based on Version 2 in Young 2013 WRR paper,
%Hypothetico-inductive data-based mechanistic modeling of hydrological 
%systems
global dailyPrecip dailyPotEvapTrans dailySS;
global t0 tEnd;

%Maximum combined contents of all stores 
cmax = parVec(1);
bexp  = parVec(2);
fQuickFlow = parVec(3);
Rs    = parVec(4);
Rq    = parVec(5);

% HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
% START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS
convFactor = 1.0;

x_loss = 0.0;

% Initialize slow tank state
%x_slow = 2.3503/(Rs*convFactor);
x_slow = 0.0;
% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0;

t= t0;
output=zeros(tEnd-t0+1,1);
counter = 1;
while t < tEnd+1
    
   Pval = dailyPrecip(t);
   PETval = dailyPotEvapTrans(t);
   
   % Compute excess precipitation and evaporation
   [UT1,UT2,x_loss] = excess(x_loss,cmax,bexp,Pval,PETval);
   
   % Partition UT1 and UT2 into quick and slow flow component
   
   UQ = fQuickFlow*(UT2 + UT1); 
   US = (1-fQuickFlow)*UT2;
   
   % Route slow flow component with single linear reservoir
   inflow = US; 
   [x_slow,outflow] = linres(x_slow,inflow,Rs); 
   QS = outflow;
   
   % Route quick flow component with linear reservoirs
   inflow = UQ; 
   k = 1; 
   while k < 4
      [x_quick(k),outflow] = linres(x_quick(k),inflow,Rq); 
      inflow = outflow; 
      k = k+1;
   end
   
   % Compute total flow for timestep
   output(counter) = (QS + outflow)*convFactor;
   t = t+1;   
   counter = counter+1;
end
sumsqr = sum((output-dailySS(t0:tEnd)).^2);
rmse = sqrt(sumsqr./(length(output)-1));
objval = rmse;