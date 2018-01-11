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
%Purpose: Generate one row in the Hermite matrix
%see also genhmat.m
function [hvec h2vec] = genhvec(N,M,P,x, polyindx)
%% 
% function genhvec(N,p)
% Input: 
%   N: the dimension of random vector
%   M: total number of terms in PC expansion
%   x: a collocation point of dimension N
%   polyindx: index of the orthogonal polynomial type
% Output: 
%    hvec, \psi_ij, 
%    h2vec, <\psi_ij^2>
%
% Note: the nth-Gaussian moments is (n-1)!
% For example <(\xi_1^2-1)(\xi_1^2-1)>
%    = <(\xi_1^4-2*\xi_1^2+1)>= 3
% Hermite Polynomials
% Order:   Form
%  H(0) =   1
%  H(1) =   x
%  H(2) =   x^2 - 1
%  H(3) =   x^3 - 3*x
%  H(4) =   x^4 - 6*x^2 + 3
% Legendre Polynomials
% Order: Form
%  H(0) =   1
%  H(1) =   x
%  H(2) =   (3*x^2 - 1)/2
%  H(3) =   (5*x^3 - 3*x)/2
%  H(4) =   (35*x^4 - 30*x^2 + 3)/8
%  H(5) =   (63*x^5 - 70*x^3 + 15*x)/8

hvec = zeros(1,M);
h2vec = zeros(1,M);
hvec(1) = 1;
h2vec(1) = 1;
% first-order terms (HP is the same as LP)
for i1=1:N
    hvec(i1+1) = getBaseFun(polyindx(i1), 1, x(i1));
    h2vec(i1+1) = getPhi2(polyindx(i1), 1);
end
% second-order terms (commented out terms are for HP)
if (P>1)
    counter=2+N;
    for i1=1:N
        for i2=i1:N
            if (i1==i2)
                hvec(counter) = getBaseFun(polyindx(i1), 2, x(i1)); % x(i1)^2-1;
                h2vec(counter) = getPhi2(polyindx(i1), 2); %2;
            else
                hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))* ...
                                getBaseFun(polyindx(i2), 1, x(i2)); %x(i1)*x(i2);
                h2vec(counter) = getPhi2(polyindx(i1), 1)* ...
                                 getPhi2(polyindx(i2), 1); %1;
            end
            counter = counter + 1;
        end
    end
end    
% third-order HP terms
if (P>2)
    for i1=1:N
        for i2=i1:N
            for i3=i2:N
                if (i1==i2 && i2==i3)
                    hvec(counter) =  getBaseFun(polyindx(i1), 3, x(i1));%x(i1)^3-3*x(i1);
                    h2vec(counter) = getPhi2(polyindx(i1), 3); %6;
                elseif (i2==i3)
                    hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))*...
                                    getBaseFun(polyindx(i2), 2, x(i2)); %x(i1)*(x(i2)^2-1);
                    h2vec(counter) = getPhi2(polyindx(i1), 1)*...
                                     getPhi2(polyindx(i2), 2); %2;
                elseif (i1==i3)
                    hvec(counter) = getBaseFun(polyindx(i2), 1, x(i2))*...
                                    getBaseFun(polyindx(i1), 2, x(i1));%x(i2)*(x(i1)^2-1);
                    h2vec(counter) =getPhi2(polyindx(i2), 1)*...
                                    getPhi2(polyindx(i1), 2); %2;
                elseif  (i1==i2)
                    hvec(counter) = getBaseFun(polyindx(i3), 1, x(i3))*...
                                    getBaseFun(polyindx(i1), 2, x(i1)); %x(i3)*(x(i1)^2-1);
                    h2vec(counter) = getPhi2(polyindx(i3), 1)*...
                                     getPhi2(polyindx(i1), 2); %2;
                else
                    hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))*...
                                    getBaseFun(polyindx(i2), 1, x(i2))*...
                                    getBaseFun(polyindx(i3), 1, x(i3)); %x(i1)*x(i2)*x(i3);
                    h2vec(counter) = getPhi2(polyindx(i1), 1)*...
                                     getPhi2(polyindx(i2), 1)*...
                                     getPhi2(polyindx(i3), 1); %1; 
                end
                counter= counter + 1;
            end
        end
    end
end
% fourth-order HP terms
if (P>3)
    for i1=1:N
        for i2=i1:N
            for i3=i2:N
                for i4=i3:N
                    if (i1==i2 && i2==i3 && i3==i4)
                        hvec(counter) = getBaseFun(polyindx(i1), 4, x(i1)); %x(i1)^4-6*x(i1)^2+3;
                        h2vec(counter) = getPhi2(polyindx(i1), 4); %24;
                    elseif (i2==i3 && i3==i4)
                        hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))*...
                                        getBaseFun(polyindx(i2), 3, x(i2)); %x(i1)*(x(i2)^3-3*x(i2));
                        h2vec(counter) = getPhi2(polyindx(i1), 1)*...
                                         getPhi2(polyindx(i2), 3); %6;
                    elseif (i1==i3 && i3==i4)
                        hvec(counter) = getBaseFun(polyindx(i2), 1, x(i2))*...
                                        getBaseFun(polyindx(i1), 3, x(i1)); %x(i2)*(x(i1)^3-3*x(i1));
                        h2vec(counter) = getPhi2(polyindx(i2), 1)*...
                                         getPhi2(polyindx(i1), 3);%6;
                    elseif (i1==i2 && i2==i4)
                        hvec(counter) = getBaseFun(polyindx(i3), 1, x(i3))*...
                                        getBaseFun(polyindx(i1), 3, x(i1)); %x(i3)*(x(i1)^3-3*x(i1));                    
                        h2vec(counter) = getPhi2(polyindx(i3), 1)*...
                                         getPhi2(polyindx(i1), 3); %6;
                    elseif (i1==i2 && i2==i3)
                        hvec(counter) = getBaseFun(polyindx(i4), 1, x(i4))*...
                                        getBaseFun(polyindx(i1), 3, x(i1)); %x(i4)*(x(i1)^3-3*x(i1));                                        
                        h2vec(counter) = getPhi2(polyindx(i4), 1)*...
                                         getPhi2(polyindx(i1), 3); %6;
                    elseif (i3==i4)
                        hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))*...
                                        getBaseFun(polyindx(i2), 1, x(i2))*...
                                        getBaseFun(polyindx(i3), 2, x(i3)); %x(i1)*x(i2)*(x(i3)^2-1);                                        
                        h2vec(counter) = getPhi2(polyindx(i1), 1)*...
                                         getPhi2(polyindx(i2), 1)*...
                                         getPhi2(polyindx(i3), 2); %2;
                    elseif (i2==i4)
                        hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))*...
                                        getBaseFun(polyindx(i3), 1, x(i3))*...
                                        getBaseFun(polyindx(i2), 2, x(i2)); %x(i1)*x(i3)*(x(i2)^2-1);                                        
                        h2vec(counter) = getPhi2(polyindx(i1), 1)*...
                                         getPhi2(polyindx(i3), 1)*...
                                         getPhi2(polyindx(i2), 2); %2;
                    elseif (i1==i4)
                        hvec(counter) = getBaseFun(polyindx(i2), 1, x(i2))*...
                                        getBaseFun(polyindx(i3), 1, x(i3))*...
                                        getBaseFun(polyindx(i1), 2, x(i1)); %x(i2)*x(i3)*(x(i1)^2-1);
                        h2vec(counter) = getPhi2(polyindx(i2), 1)*...
                                         getPhi2(polyindx(i3), 1)*...
                                         getPhi2(polyindx(i1), 2); %2;
                    elseif (i1==i3)
                        hvec(counter) = getBaseFun(polyindx(i2), 1, x(i2))*...
                                        getBaseFun(polyindx(i4), 1, x(i4))*...
                                        getBaseFun(polyindx(i1), 2, x(i1)); %x(i2)*x(i4)*(x(i1)^2-1);
                        h2vec(counter) = getPhi2(polyindx(i2), 1)*...
                                         getPhi2(polyindx(i4), 1)*...
                                         getPhi2(polyindx(i1), 2); %2;
                    elseif (i2==i3)
                        hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))*...
                                        getBaseFun(polyindx(i4), 1, x(i4))*...
                                        getBaseFun(polyindx(i2), 2, x(i2)); %x(i1)*x(i4)*(x(i2)^2-1);
                        
                        h2vec(counter) = getPhi2(polyindx(i1), 1)*...
                                         getPhi2(polyindx(i4), 1)*...
                                         getPhi2(polyindx(i2), 2); %2;
                    elseif (i1==i2)                    
                        hvec(counter) = getBaseFun(polyindx(i3), 1, x(i3))*...
                                        getBaseFun(polyindx(i4), 1, x(i4))*...
                                        getBaseFun(polyindx(i1), 2, x(i1));%x(i3)*x(i4)*(x(i1)^2-1);
                        h2vec(counter) = getPhi2(polyindx(i3), 1)*...
                                         getPhi2(polyindx(i4), 1)*...
                                         getPhi2(polyindx(i1), 2); %2;
                    else
                        hvec(counter) = getBaseFun(polyindx(i1), 1, x(i1))*...
                                        getBaseFun(polyindx(i2), 1, x(i2))*...
                                        getBaseFun(polyindx(i3), 1, x(i3))*...
                                        getBaseFun(polyindx(i4), 1, x(i4)); %x(i1)*x(i2)*x(i3)*x(i4);
                        h2vec(counter) = getPhi2(polyindx(i1), 1)*...
                                         getPhi2(polyindx(i2), 1)*...
                                         getPhi2(polyindx(i3), 1)*...
                                         getPhi2(polyindx(i4), 1); %1;
                    end
                    counter= counter + 1;
                end
            end
        end
    end
end
end
function [res] = getBaseFun(itype, porder, x)
%  H(0) =   1
%  H(1) =   x
%  H(2) =   x^2 - 1
%  H(3) =   x^3 - 3*x
%  H(4) =   x^4 - 6*x^2 + 3
% Legendre Polynomials
% Order: Form
%  H(0) =   1
%  H(1) =   x
%  H(2) =   (3*x^2 - 1)/2
%  H(3) =   (5*x^3 - 3*x)/2
%  H(4) =   (35*x^4 - 30*x^2 + 3)/8
%  H(5) =   (63*x^5 - 70*x^3 + 15*x)/8
    if (itype==1)
        switch porder
            case 0
                res = 1;
            case 1
                res = x;
            case 2
                res = x^2-1;
            case 3
                res = x^3-3*x;
            case 4
                res = x^4-6*x^2+3;
        end
    elseif (itype==2)
        switch porder
            case 0
                res = 1;
            case 1
                res = x;
            case 2
                res = 0.5*(3*x^2-1);
            case 3
                res = 0.5*(5*x^3-3*x);
            case 4
                res = 0.125*(35*x^4-30*x^2+3);
        end
    end
end
function [res] = getPhi2(itype, porder)
    if (itype==1)
        switch porder
            case 1
                res = 1;
            case 2
                res = 2;
            case 3
                res = 6;
            case 4
                res = 24;
        end
    elseif (itype==2)
        switch porder
            case 1
                res = 1/3;
            case 2
                res = 1/5;
            case 3
                res = 1/7;
            case 4
                res = 1/9;
        end
    end
end