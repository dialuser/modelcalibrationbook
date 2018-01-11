function x = rtlsqep(A,b,L,delta,x0,typeL)
% RTLSQEP solves the regularized TLS problem
%             || Ax - b ||^2
%        min ---------------, subject to ||Lx||^2 = delta
%         x    ||x||^2 + 1
% 
% Input args:
%    A - mxn matrix of coefficients
%    b - mx1 right-hand-side vector 
%    L - nxn regularization matrix, full-rank
%    delta - regularization constant, i.e., 
%            the constraint ||Lx||^2 = delta is imposed
%    x0 - initial guess for the solution (can be empty)
%    typeL - specifies special structure of L:
%         'id' - identity matrix; 'der1' - first oder derivative; 
%         'der2' - second order derivative; 'full' - full matrix (default) 
%
% Output args:
%    x - solution of the RTLS problem
%
% Method:
% Iteratively, a quadratic eigenvalue problem is solved,
% equivalent to solving a system of the form:
%
%  (B(x_old) + lam Q)x = d(x_old),     x'Qx = delta,  where Q = L'*L.
%
%  Only largest quadratic eigenvalue and corresponding eigenvector
%  are needed; they are computed with Matlab's eigs function (from ARPACK).

% Contributor: Diana Sima, KU Leuven, october 2003.

  
% dimensions and arguments check
if nargin<4, error('There should be at least 4 input arguments'), end

[m,n] = size(A); 

[mb,nb] = size(b);
if (nb~=1) 
  error('Only column vector RHS b is accepted.')
end
if (mb~=m) 
  error('b should have as many elements as number of rows in A')
end

[mL,nL] = size(L);
if (mL~=nL) | (nL~=n) 
  error('L should be a square matrix of size n, the column dimension of A')
end

if ~isscalar(delta) | (delta<=0)
  error('delta should be a positive scalar')
end

if nargin<5 | isempty(x0)
  x0 = rand(n,1); 
elseif (size(x0,1)~=n) | (size(x0,2)~=1)
  error(['x0 should be a column vector of length ',num2str(n)])
end
if nargin<6 | isempty(strfind(' id der1 der2 full ',[' ',typeL,' ']))
  typeL = 'full'; 
end

% set tolerance level
tol = 1e-8;

% set maximum number of iterations
maxiter = 10; 

% initializations
x  = x0; 
mu = x'*x + 1;
ro = norm(A*x-b)^2/mu;
h  = (b'*A)'; h = solveL(L,h,typeL,'trans'); 
opts.disp = 0; opts.tol = tol; opts.maxit = 5; opts.p = 25; 
err = 1; count = 1; 

% main loop
while (err>tol & count<maxiter) 

  count = count + 1;
  x1    = x;
   
  % solve eigenvalue problem
  [u,l,flag] = eigs(@computeB,2*n,1,'LR',opts,A,b,L,h/mu,ro,mu, ...
                    delta,typeL);  
  if flag~=0, warning('eigenvalue problem not solved accurately'), end 
  
  % compute solution from eigenvector
  opts.v0 = u;
  u(1:n)  = []; 
  u = u*delta*mu/(h'*u); 
  x = l*u + computeW(u,A,L,ro,mu,typeL); 
  x = solveL(L,x,typeL); 
  mu  = x'*x + 1;
  ro  = norm(A*x-b)^2/mu;
  err = norm(x1-x)/norm(x1);

end

% end rtlsqep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auxiliary functions: %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = computeB(x,A,b,L,h,ro,mu,delta,typeL)
% compute matrix vector product y = B*x, where
% B = [-2W, -W*W + hh'/(del^2);
%        I, 0]  
% W*x is computed with 'computeW'
[m,n] = size(A);
v = x(1:n); u = x(n+1:2*n); y2 = v;
y1 = computeW(v,A,L,ro,mu,typeL); y1 = -2*y1; 
z = computeW(u,A,L,ro,mu,typeL);
z = computeW(z,A,L,ro,mu,typeL);
z = -z + (h'*u)/delta*h; 
y1 = y1 + z; 
y = [y1;y2]; 
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Wx = computeW(x,A,L,ro,mu,typeL)
% matrix vector product:
% Wx := W*x, where
%                W = (L^(-T)A'AL^(-1) - ro (L'L)^(-1))/mu

z1 = solveL(L,x,typeL);         % z1 = L\x;
z1 = A*z1; z1 = (z1'*A)';       % z1 = A'Az1;
z1 = solveL(L,z1,typeL,'trans');% z1 = L'\z1;
%
z2 = solveL(L,x,typeL);         % z2 = L'\x;
z2 = solveL(L,z2,typeL,'trans');% z2 = L\z2;
z2 = -ro*z2;
%
Wx = (z1+z2)/mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = solveL(L,x,typeL,trans)
% solve the system Lz = x or L'z = x
% using fast methods if L has special structure
% 
% typeL = 'id' | 'der1' | 'der2' |  'full'  

n = length(x); z = zeros(n,1);
if nargin<4, trans ='notrans'; end
if nargin<3, typeL ='full'; end
if strcmp(trans,'notrans')==1
  switch typeL 
    case 'id'
        z = x;
    case 'der1'
        z(1) = x(1);
        for j = 2:n
            z(j) = z(j-1) + x(j);
        end
    case 'der2'
        z(1) = x(1); z(2) = 2*z(1) + x(1);
        for j = 3:n
            z(j) = 2*z(j-1) - z(j-2) +  x(j);
        end
    case 'full'
        z = (x'/L)';
  end
else
   switch typeL 
    case 'id'
        z = x;
    case 'der1'
        z(n) = x(n);
        for j = n-1:-1:1
            z(j) = z(j+1) + x(j);
        end
    case 'der2'
        z(n) = x(n); z(n-1) = 2*z(n) + x(n-1);
        for j = n-2:-1:1
            z(j) = 2*z(j+1) - z(j+2) +  x(j);
        end
    case 'full'
        z = L\x;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%