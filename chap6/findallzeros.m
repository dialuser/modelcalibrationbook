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
%Used by example6_5.m
%====================================================================
function  x  =  findallzeros( fhandle, xguess, maxroots)
%
%    x = findallzeros( fhandle, xguess, maxroots)
%   ==================================
%
%   attemps to find multiple roots of a function.
%
%  fhandle = function handle, prefix by @
%  xguess = a scalar value used as initial guess
%           ( = 0   if not supplied )
% maxroots = the number of roots needed
%           ( = 1   if not supplied )

%   (c) 2005 by Tomas Co
%   @ Michigan Technological University

    x   = [ ]                ;
    
    if nargin==1
        xguess =0;
        maxroots = 1;
    end
    
    if nargin==2
        maxroots  = 1;
    end

    nroots = 0;
    
    while nroots<maxroots
        if  abs(feval(fhandle,xguess))<=1e-16,
            xroot=xguess;
        else
            [xroot,val,flag] = fzero(fhandle,xguess);
            if flag <= 0
                return
            else
                if (nroots>0 && isempty(find(abs(x-xroot)<1e-8)))
                    x=[x; xroot];
                    nroots=nroots+1;
                    xguess = xroot;                    
                elseif (nroots==0 && xroot>0.1)
                    x=xroot;
                    nroots=nroots+1;
                end
            end
        xguess=xguess+0.05;
        end
    end