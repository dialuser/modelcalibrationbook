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
%%
%Author: Alex Sun
%Purpose: Generate hermite polynomial matrix for response surface
%Date: 5/5/2011
%Revised: 3/17/2012
%Revise for generalized chaos expansion
%% 
function [hmat,xmat,h2mat] = genhmat(N, P, polyindx)
%function genhmat(N,P)
%Input:
%   N: dimension of the random input
%   P: order of the polynomial to be used
%  polyindx: indicator of the polynomial
%       1, Hermite polynomial
%       2, Legendre polynomial       
%Ouput:
%   hmat: the coefficient matrix for PCM
%Reference:
% Eqs (2) and (7) in 
% S. Oladyshkin, H. Class, R Helmig, and W. Nowak (2011)
% An integrative approach to robust design and probabilistic risk 
% assessment for CO2 storage. Comput Geosci  
% Li and Zhang, WRR 2007 paper
%  PCM for flow in porous media
%% M is the number of terms in PC expansion, see Eq (2)
addpath('c:/alex/cO2/orthopoly');
P1=P+1;
M = factorial(N+P)/(factorial(N)*factorial(P));
hmat = zeros(M,M);
h2mat = hmat;
%% Generating collocation points of using roots of (P+1) HP
%get Hermite roots
r = HermiteRoots(P+1);
%get Legendre roots
rL = LegendreRoots(P+1);
%sort roots according to probability
r = sortme(r);
rL = sortme(rL);

%Initialize the collation point matrix, xmat w/ the first root in root vect
%which has the biggest probability 
%(for uniform, I assume larger number corresponds to higher probabiliy
% smaller number for smaller probability; therefore, the order should be
% biggest, smallest, bigger, smaller, ....
xmat = zeros(M,N);
for i=1:N
    if (polyindx(i)==1)
        xmat(:, i) = r(1);
    elseif (polyindx(i)==2)
        xmat(:, i) = rL(1);
    end
end
        
%Algorithm: 
% the first row consists of the smallest collocation point r(1)
% starting from the second row, fill the matrix with r(2), first with one
% element and then two elements, and so on. The width is controlled by iw
% 
[hmat(1,:) h2mat(1,:)]= genhvec(N,M,P,xmat(1,:),polyindx);
counter = 2;
rvec0 = zeros(1,N);
for i=1:N
    if (polyindx(i)==1)
        rvec0(i) = rvec0(i)+r(1);
    elseif (polyindx(i)==2)
        rvec0(i) = rvec0(i)+rL(1);
    end
end
for k=2:P1
    for i=1:N
        rvec = rvec0;
        rvec(i) = getRoot(i,k); %(r(k); 
        [hmat h2mat xmat counter]=checkrank(N, M, P, rvec, hmat, h2mat, xmat, counter, polyindx);   
        if (counter>M) 
            return;
        end
    end
    for i1=1:N-1 
        for k1=2:P1
        for i2=i1+1:N
            rvec=rvec0;
            rvec(i1)=getRoot(i1,k);%r(k);
            rvec(i2)=getRoot(i2,k1);%r(k1);  
            [hmat h2mat xmat counter]=checkrank(N, M, P, rvec, hmat, h2mat, xmat, counter, polyindx);
            if (counter>M) 
                return;
            end        
        end
        end
    end
    for i1=1:N-2
        for k1=2:P1
        for i2=i1+1:N-1
            for k2=2:P1
            for i3=i2+1:N
                rvec=rvec0;
                rvec(i1)=getRoot(i1,k);%r(k);
                rvec(i2)=getRoot(i2,k1); %r(k1);
                rvec(i3)=getRoot(i3,k2); %r(k2);
                [hmat h2mat xmat counter]=checkrank(N, M, P, rvec, hmat, h2mat, xmat, counter, polyindx);
                if (counter>M) 
                    return;
                end                        
            end
            end
        end
        end
    end
    for i1=1:N-3
        for k1=2:P1
        for i2=i1+1:N-2
            for k2=2:P1
            for i3=i2+1:N-1
                for k3=2:P1
                for i4=i3+1:N
                    rvec=rvec0;
                    rvec(i1)=getRoot(i1,k);%r(k);
                    rvec(i2)=getRoot(i2,k1);%r(k1);
                    rvec(i3)=getRoot(i3,k2); %r(k2);
                    rvec(i4)=getRoot(i4,k3); %r(k3);
                    [hmat h2mat xmat counter]=checkrank(N, M, P, rvec, hmat, h2mat, xmat, counter,polyindx);
                    if (counter>M) 
                        return;
                    end                            
                end
                end
            end
            end
        end
        end
    end
end 
    function  [val] = getRoot(ii, kk)
        if (polyindx(ii)==1)
            val = r(kk);
        elseif(polyindx(ii)==2)
            val = rL(kk);
        end
    end
end
function [hmat h2mat xmat counter]=checkrank(N, M, P, rvec, hmat, h2mat, xmat, counter, polyindx)
    [hvec, h2vec] = genhvec(N,M,P,rvec, polyindx);
    tempmat = [hmat(1:counter-1,:); hvec];
    if (rank(tempmat)==counter )
        hmat(counter, :) = hvec;
        h2mat(counter,:) = h2vec;
        xmat(counter, :) = rvec;
        disp(counter);
        counter = counter + 1;                
    end
end
function [res] = sortme(r)
    d = length(r);
    res =zeros(d,1);
    if (d/2==fix(d/2))
        mid = d/2;
        res(1) = r(mid);
        counter = 2;
        for jj=1:mid
            res(counter) = r(mid+jj);
            counter=counter+1;
            if (counter>=d) 
                break; 
            end
            res(counter) = r(mid-jj);
            counter=counter+1;
        end
    else
        mid = ceil(d/2);
        res(1) = r(mid);
        counter=2;
        for jj=1:mid-1
            res(counter)=r(mid-jj);
            res(counter+1)=r(mid+jj);
            counter=counter+2;
            if (counter>=d)
                break;
            end
        end
    end
end
