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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test1D
%Case 3
%Author: Alex Sun
%Date: 3/23/2005
%Project:Source ID
%Generates a synthetic plume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obsvals]=samplegen1d3(prdconcs)
%Params:
%prdconcs, source concentration during each release prd
%
%Global variables
global obswells obstimes nwells ntimes nsamples;
%
obsvals=[];
for j=1:ntimes
    %call plumeGen1d3 to generate concentration plume
    [concres]=plumegen1d3(obstimes(j), prdconcs);
    totalnodes=size(concres,1);
    %sample the plume
    for i=1:totalnodes
        for k=1:nwells
            if(obswells(k) == concres(i,1))
                obsvals=[obsvals concres(i,2)];
            end
        end
    end
end