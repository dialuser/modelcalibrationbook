%adapted from
%DREAM and 
%https://github.com/NLeSC/esibayes/tree/master/examples/scemua-so-hymod-batch/model

function [x_f, outflow] = linres(x_f,inflow,R)
% Linear reservoir
x_f = (1-R)*x_f + (1-R)*inflow;
outflow = (R/(1-R))*x_f;

