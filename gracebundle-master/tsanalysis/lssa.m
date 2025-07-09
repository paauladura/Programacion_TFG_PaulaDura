function [c,s,m,t] = lssa(x,y,f,dy,mn,trnd)

% LSSA computes the least-squares spectrum of a given dataset for the
% prescribed frequencies.
%
% F = lssa(x,y,f)
% F = lssa(x,y,f,mn,trnd)
%
% INPUT
% x     -   observation points
% y     -   observations
% dy    -   precision of observations. either standard deviations or the 
%           full variance-covariance matrix.
% f     -   frequencies whose coefficients have to be estimated, and these
%           frequencies must be given as non-anugular frequencies.
% mn    -   switch that toggles the estimation of a mean: '0' [default] for
%           no mean and '1' otherwise.
% trnd  -   toggle switch for trend estimation: '0' [default] and '1'
%
% OUTPUT
% [c s] -   Vector of cosine and sine coefficients of the frequencies
% m     -   mean
% t     -   trend
%--------------------------------------------------------------------------
% See also lombscargle fft
%--------------------------------------------------------------------------

if (min(size(x)) ~= 1) || (min(size(y)) ~= 1) || (min(size(f)) ~= 1)
    error('x,y and f must all be vectors');
end

x = x(:);
y = y(:);
f = f(:);

if nargin < 3
    error('Insufficient input arguments')
elseif nargin == 3
    dy = ones(size(y));
    mn   = false;
    trnd = false;
elseif nargin == 4
    mn   = false;
    trnd = false;
elseif nargin == 5
    mn   = logical(mn);
    trnd = false;
elseif nargin == 6
    mn   = logical(mn);
    trnd = logical(trnd);
end



[r,c] = size(dy);
if max(r,c) ~= max(size(y))
    error('The observations and their precision must be of the same length')
elseif (min(r,c) ~= 1) | (r ~= c)
    error('Please provide either the precision of the observations or its covariance matrix')
elseif min(r,c) == 1
    P = diag(1./(dy.^2));
elseif r == c
    P = dy\eye(size(dy));
end


flen = length(f);

arg  = 2*pi*(x)*f(:)';

A = [cos(arg) sin(arg)];

if trnd
    A = [x A]; 
else
    t = NaN;
end
if mn
    A = [ones(size(x)) A]; 
else
    m = NaN;
end

N = A' * P;
y = N * y;
N = N * A;

Q = N\eye(size(N));

F = Q * y;
if mn, m = F(1); F = F(2:end); end
if trnd, t = F(1); F = F(2:end); end
c = F(1:flen);
s = F((1+flen):end);

tmp = Q;

if mn && trnd
    Q.mm = tmp(1,1);
    Q.mt = tmp(2,1);
    Q.mc = tmp(3:flen+2,1);
    Q.ms = tmp(flen+3:end,1);

    Q.tt = tmp(2,2);
    Q.tc = tmp(3:flen+2,2);
    Q.ts = tmp(flen+3:end,2);

    tmp  = tmp(3:end,3:end);
elseif mn && ~trnd
    Q.mm = tmp(1,1);
    Q.mc = tmp(2:flen+1,1);
    Q.ms = tmp(flen+2:end,1);

    tmp  = tmp(2:end,2:end);
elseif ~mn && trnd
    Q.tt = tmp(1,1);
    Q.tc = tmp(2:flen+1,1);
    Q.ts = tmp(flen+2:end,1);

    tmp = tmp(2:end,2:end);
end

Q.cc = tmp(1:flen,1:flen);
Q.sc = tmp(flen+1:end,1:flen);
Q.ss = tmp(flen+1:end,flen+1:end);
