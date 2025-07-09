function ip = iplmquad(l, m, theRAD, dt)

% Integrals of the fully normalized associated Legendre functions computed 
% with adaptive Lobatto quadrature. Matlab function QUADL implements a high 
% order method using an adaptive Gauss/Lobatto qudrature rule for 
% integration. l & m are degree and order of the respective Plm function and
% theRAD is the theta-vector given in radians. The integration span dt is an
% optional argument (default: theRAD(i+1) - theRAD(i))
%
% HOW:     IP  = IPLMQUAD(l, m, theRAD, dt)
%          IP  = IPLMQUAD(l, m, theRAD)    (dt = theRAD(i+1) - theRAD(i))
%
% USES:
%    plm

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Dimitris TSOULIS (DT), IAPG, TU-Munich  
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-02-13: MR, input argument in radians, 
%                    replacing obsolete 'quad8' by 'quadl',
%                    replacing function name by function handle and putting 
%                    integral kernel function as internal function into this 
%                    file, updated help text.
%    2012-01-23: MA, intpnm_internal: input of plm in radian
%    1999-01-??: DT, initial version
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

if ~exist('dt', 'var')
    dt = theRAD(2) - theRAD(1);
end

for t = 1:length(theRAD)
    ip(t) = quadl(@intpnm_internal, theRAD(t) - dt/2, theRAD(t) + dt/2, [], [], l, m);
end

%% -------------------------------------------
% INTPNM_INTERNAL Integral kernel of the integrals of the fully normalized 
%        associated Legendre functions.
function i = intpnm_internal(theRAD, n, m)
i = plm(n, m, theRAD) .* sin(theRAD);


