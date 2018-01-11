%adapted from
%DREAM and
%https://github.com/NLeSC/esibayes/tree/master/examples/scemua-so-hymod-batch/model
%https://github.com/NLeSC/esibayes/blob/master/manual/other/excess.m

function [OV,xHuz] = excess(xHuz,Cpar,bexp,Pval,PETval,Huz, Kv)
% this function calculates excess precipitation and evaporation
%
%Inputs:
%     Huz, Maximum height of soil moisture accounting tank - Range [0, Inf]
%     bexp, Scaled distribution function shape parameter   - Range [0, 2]
%     fQuickFlow, Quick/slow split parameter               - Range [0, 1]
%     Rq, Quickflow routing tanks' rate parameter          - Range [0, 1]
%     Rs, Slowflow routing tank's rate parameter           - Range [0, 1]
%     Kv, Vegetation adjustment to PE                      - Range [0, 2]
%Outputs:
%     OV1, excess rainfall 
%     OV2, excess 
%xHuz, current content

%Storage contents at begining
Cbeg = Cpar*(1.0 - power(1.0-(xHuz/Huz),1.0+bexp));
%Compute effective rainfall filling all storage elements
OV2 = max(0.0, Pval + xHuz - Huz);
%Remaining net rainfall
PPinf = Pval - OV2;
%New actual height
Hint = min(Huz, xHuz+PPinf);
%New storage content
Cint = Cpar*(1.0-power(1.0-(Hint/Huz),1.0+bexp));
%Additional effective rainfall produced by stores smaller than Cmax
OV1 = max(0.0, PPinf + Cbeg - Cint);
%Compute total effective rainfall
OV = OV1 + OV2;
%Computer actual evapotranspiration
AE = min(Cint,(Cint/Cpar)*PETval*Kv);
%Storage contents after ET
xCuz = max(0.0,Cint - AE);
%Compute final height of the reservoir
xHuz = Huz*(1.0-power(1.0-(xCuz/Cpar),1.0/(1.0+bexp)));
