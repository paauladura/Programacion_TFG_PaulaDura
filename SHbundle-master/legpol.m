function  p = legpol(l, m, thetaRAD)

% LEGPOL creates a matrix of unnormalized Legendre functions 
%
% HOW: p = legpol(l,thetaRAD)
%      p = legpol(l,m,thetaRAD)
%
% IN:
%    l ........ vector of degree values. Integer, but not necessarily monotonic.
%               For l < m a vector of zeros will be returned.
%    m ........ selected order (scalar). If absent, m=0 is assumed.
%    thetaRAD ... vector of colatitudes [rad]
%
% OUT:
%    p ........ matrix of unnormalized Legendre functions with length(thetaRAD) rows 
%               and length(l) columns. In case l or thetaRAD is scalar, the output vector
%               will follow the shape of respectively thetaRAD or l. 
% 
% SEE ALSO:
%    PLM, YLM
%
% REMARKS:
%    For m>0, we have Legendre functions, not polynomials. The m-file name is 
%    therefore wrong.

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias WEIGELT (MW), DoGE, UofC
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-01-29: MA, comments
%    2013-01-23: MW, input argument thetaRAD [deg -> radian]
%    1999-03-17: NS, help text brushed up
%    1994-08-17: NS, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

% Some input checking.
if nargin == 2
   thetaRAD = m;
   m = 0;
elseif nargin ~= 3
   error('LEGPOL requires 2-3 input arguments.')
end
if min(size(l)) > 1,   error('Degree L must be vector (or scalar)'), end
if numel(m) > 1,  error('Order M must be scalar.'), end
if sum(rem(l,1)) ~= 0, error('Degree L must be integer'), end
if rem(m,1) ~= 0,      error('Order M must be integer'), end

% Preliminaries.
[lrow,lcol] = size(l);
[trow,tcol] = size(thetaRAD);
lmax = max(l);
if lmax < m, error('Largest degree still smaller than order m.'), end
n    = length(thetaRAD);			% number of latitudes
t    = thetaRAD(:);			    % column vector in [rad]
x    = cos(t);
y    = sin(t);
mfix = m;					% m can be used now as running index.
lvec = l(:)';					% l can be used now as running index.
if min(t)<0 || max(t)>pi
    warning('Is the co-latitude ''thetaRAD'' given in radian?')
end

% Recursive computation of the temporary matrix ptmp, containing the Legendre
% functions in its columns, with progressing degree l. The last column of
% ptmp will contain zeros, which is useful for assignments when l < m.
ptmp = zeros(n,lmax-mfix+2);

% sectorial recursion
p = ones(n,1);
for m = 1:mfix
   p = (2*m-1) * p .* y;
end
ptmp(:,1) = p;				% The 1st column of ptmp.

% l-recursion
pold = zeros(n,1);
for l = mfix+1:lmax
   col    = l - mfix + 1;			% points to the next column of ptmp
   pnew   = (2*l-1)/(l-mfix) * x .* p - (l+mfix-1)/(l-mfix) * pold;
   pold   = p;
   p      = pnew;
   ptmp(:,col) = p;
end


% The Legendre functions have been computed. What remains to be done, is to
% extract the proper columns from ptmp, corresponding to the vector lvec. 
% If l or thetaRAD is scalar the output matrix p reduces to a vector. It should
% have the shape of respectively thetaRAD or l in that case.
p    = zeros(n, length(lvec));		% size declaration.
lind  = find(lvec < mfix);			% index into l < m
pcol  = lvec - mfix + 1;			% index into columns of ptmp
pcol(lind) = (lmax-mfix+2)*ones(size(lind));	% Now l < m points to last col.
p    = ptmp(:,pcol);			% proper column extraction 


% reshape output P in case L or THETARAD are scalar
if numel(lvec) == 1,  p = reshape(p,trow,tcol); end
if numel(thetaRAD) == 1, p = reshape(p,lrow,lcol); end
