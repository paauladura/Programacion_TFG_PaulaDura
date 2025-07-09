function clm = klm2clm(klm)

% KLM2CLM converts spherical coefficients ordered in [l (-/+)m Klm] to
% [l m Clm Slm].
%
% clm = klm2clm(klm)
%
% INPUT
% klm   - Spherical harmonic coefficients given in a look-up-table, in 
%         which the coefficients are ordered
%           [0    0     C_(0,   0); 
%            1   -1     S_(1,  -1);
%            1    0     C_(1,   0);
%            1    1     C_(1,   1);
%            2   -2     S_(2,  -2);
%                ...
%            N  N-1     C_(N, N-1);
%            N    N     C_(N,   N);
%
% OUTPUT
% clm    - Spherical harmonic coefficients arranged in a look-up-table
%          and ordered in order-leading format.
%           [0      0   C_(0, 0)    S_(0, 0);
%            1      0   C_(1, 0)    S_(1, 0);
%            2      0   C_(2, 0)    S_(2, 0);
%                       ...
%            n      m   C_(n, m)    S_(n, m);
%                       ...
%            N      N   C_(N, N)    S_(N, N)];
%
%-----------------------------------------------------------------------

% Balaji Devaraju, IfE, Hannover
% Project: SHbundle

% Version control
%   BD  2016-09-30  Initial version
%-----------------------------------------------------------------------



narginchk(1, 1)

lmax = max(klm(:,1));
L   = lmax + 1;
N   = (L+1)*L/2;
M   = cumsum(0:lmax)';


clm = zeros(N, 4);

idx = (klm(:,2) < 0);

cidx = L*klm(~idx,2)     + klm(~idx, 1)  + 1 - M(klm(~idx,2)+1);
sidx = L*abs(klm(idx,2)) + klm(idx, 1)   + 1 - M(abs(klm(idx,2))+1);

clm(cidx,1:3) = klm(~idx,1:3);
clm(sidx,4) = klm(idx,3);
