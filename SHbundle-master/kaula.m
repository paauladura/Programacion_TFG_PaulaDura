function degrms = kaula(degree,sf)

% KAULA(L,SF) computes Kaula's rule of thumb. 
%
% Kaula's rule represents the expected size (rms-sense) of a single
% spherical harmonic coefficient.
%
% IN:
%    L(degree) ..... can be scalar, vector or matrix.
%    SF(sf) ........ scale factor         (Default is chosen as sqrt(0.5))
%
% USES:
%    uberall/isint, uberall/replace
%
% SEE ALSO:
%    TSCHERRAPP

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico Sneeuw (NS), IAPG, TU-Munich
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-02-13: MR, change function names, brush up comments
%    1998-06-03: NS, Inf output in case of degree=0
%                    no checking for negative degree anymore
%    1994-04-13: NS, initial version
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

if ~isint(degree), error('The degree vector contains non-integers.'), end

if nargin == 1
   sf = sqrt(0.5);
end
degrms = sf * 1e-5 * degree.^(-2);

degrms = replace(degrms, degree < 2, Inf);
