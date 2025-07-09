function klm = clm2klm(clm)

% CLM2KLM converts spherical coefficients ordered in [l m Clm Slm] to
% [l (-/+)m Klm].
%
% klm = clm2klm(klm)
%
% INPUT
% clm    - Spherical harmonic coefficients arranged in a look-up-table
%          and ordered in the following format
%           [0      0   C_(0, 0)    S_(0, 0);
%            1      0   C_(1, 0)    S_(1, 0);
%            2      0   C_(2, 0)    S_(2, 0);
%                       ...
%            n      m   C_(n, m)    S_(n, m);
%                       ...
%            N      N   C_(N, N)    S_(N, N)];
%          The coefficients need not be ordered in order-leading format.
%
% OUTPUT
% klm   - Spherical harmonic coefficients given in a look-up-table, in 
%         which the coefficients are ordered as
%           [0    0     C_(0,   0); 
%            1   -1     S_(1,  -1);
%            1    0     C_(1,   0);
%            1    1     C_(1,   1);
%            2   -2     S_(2,  -2);
%                ...
%            N  N-1     C_(N, N-1);
%            N    N     C_(N,   N);
%
%-----------------------------------------------------------------------

% Balaji Devaraju, IfE, Hannover
% Project: SHbundle

% Version control
%   BD  2017-06-26  Initial version
%-----------------------------------------------------------------------



narginchk(1, 1)

lmax = max(clm(:,1));
L   = lmax + 1;
N   = L^2;
M   = cumsum(0:lmax)';

m   = clm(:,2);
l   = clm(:,1);
idx = (m ~= 0);
cid = l.^2      + (l      + 1) + m;
sid = l(idx).^2 + (l(idx) + 1) - m(idx);


klm = zeros(N, 3);
klm(cid, 1) = l;
klm(sid, 1) = l(idx);

klm(cid, 2) = m;
klm(sid, 2) = -m(idx);

klm(cid, 3) = clm(:, 3);
klm(sid, 3) = clm(idx, 4);
