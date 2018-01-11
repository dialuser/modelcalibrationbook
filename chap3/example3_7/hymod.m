%adapted from
%DREAM and
%https://github.com/NLeSC/esibayes/tree/master/examples/scemua-so-hymod-batch/model
function [objval output] = hymod(parVec)
global dailyPrecip dailyPotEvapTrans dailySS;
global t0 tEnd;

%Maximum combined contents of all stores 
Huz = parVec(1);
bexp  = parVec(2);
fQuickFlow = parVec(3);
Rs    = parVec(4);
Rq    = parVec(5);
Kv    = 1.0;

% HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
% START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS

% Initial height corresponding to SMA tank contents
x_Huz = 0.0;
% Initialize slow tank state
%x_slow = 2.3503/(Rs*convFactor);
x_slow = 0.0;
% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0;

cPar = Huz/(1.0+bexp);

t= t0;
output=zeros(tEnd-t0+1,1);
counter = 1;
while t < tEnd+1
    
   Pval = dailyPrecip(t);
   PETval = dailyPotEvapTrans(t);
   
   % Compute excess precipitation and evaporation
   [UT,x_Huz] = excess(x_Huz,cPar,bexp,Pval,PETval,Huz,Kv);
   
   % Partition UT1 and UT2 into quick and slow flow component
   UQ = fQuickFlow*UT; 
   US = (1-fQuickFlow)*UT;
   
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
   output(counter) = QS + outflow;
   t = t+1;   
   counter = counter+1;
end
sumsqr = sum((output-dailySS(t0:tEnd)).^2);
meanobs = mean(dailySS(t0:tEnd));
rmse = sqrt(sumsqr./(length(output)-1));
%nsc = -(1- sumsqr/sum((dailySS(t0:tEnd)-meanobs).^2));
nsc = sumsqr/sum((dailySS(t0:tEnd)-meanobs).^2);
objval = rmse;