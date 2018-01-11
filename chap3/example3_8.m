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
%Date: $20120102
%Example 3.8
%Purpose: Use the USGS/LOADEST Example_1 data to illustrate GA
%Dependency:
%	Requires Matlab Global Optimization Toolbox
%====================================================================
function [] = example3_8()
%Read data from /data folder
fid = fopen('data/calib.inp', 'r');
%skip the header
for i=1:16
    fgets(fid);
end
res = textscan(fid, '%s %*n %n %n');
cdatesstr = res{1}(:);
flowrate = res{2}(:); %[ft^3/s]
conc = res{3}(:);     %[mg/l], 1L = 0.035315 ft^3
cdates = datenum(cdatesstr, 'yyyymmdd');
fclose(fid);

figure(1), clf;
%Need to convert from conc to kg/ft^3 and Q to ft^3/d
obsload = log(conc.*flowrate*86400*1e-6/0.035315); %[kg/d]

[AX,H1,H2] = plotyy(cdates,obsload, cdates, flowrate);

set(get(AX(1),'Ylabel'),'String','Phosphorus Load (kg/d)') 
set(get(AX(2),'Ylabel'),'String','Flow rate (ft^3/d)') 
set(AX(1), 'XTick', cdates, 'YColor', 'k', 'YScale', 'log');
set(AX(2), 'XTick', cdates, 'YColor', 'k', 'YScale', 'log', 'YLim', [1000 1e5]);

set(H1,'Marker','o', 'Color', 'r');
set(H2,'Marker','x', 'Color', 'b');

datetick(AX(1),'x', 'yyyy', 'keeplimits');
datetick(AX(2),'x', 'yyyy', 'keeplimits');
legend('Load', 'Flowrate');
options = gaoptimset('MutationFcn',@mutationadaptfeasible, ...
                    'Generations', 100, ...
                    'Display', 'diagnose', ...
                    'TolFun', 1e-7, ...
                    'CrossoverFraction', 0.8, ...
                    'PlotFcns', @gaplotbestf)
 

lb = [0; 0];
[thetaparam,fval,exitflag,output] = ga(@loadmodel,2,[],[],[],[],...
lb,[],[],options);

disp('Successfully finished Example 3.8');
    function [res] = loadmodel(alpha)
        %calculate log load
        %Get centered streamflow (see section 2.3 in loadest manual)
        lnQ = log(flowrate);
        lnQm = mean(lnQ);
        lnQc = lnQm+sum((lnQ-lnQm).^3)./(2*sum((lnQ-lnQm).^2));
        lnQ = lnQ-lnQc;
        
        estload = alpha(1)+alpha(2)*lnQ; %[kg/d]
        
        res = norm(estload-obsload);
    end
end
