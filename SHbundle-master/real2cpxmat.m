function C = real2cpxmat(n)

% REAL2CPXMAT generates a real to complex conversion of spherical harmonics
% for the given degree.
%
% IN:
%    n ....... Degree of spherical harmonics (integer)
% 
% OUT:
%    C ....... Conversion matrix for the given degree.
%
% SEE ALSO:
%    real2cpxsh

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2010-12-09: BD, initial version
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

m = (1:n)';

C = diag([ones(n,1)*1i; 1/sqrt(2); (-1).^(m)]);
C = C + fliplr(diag([ones(n,1);1/sqrt(2); (-1).^(m)*(-1i)]));
C = C/sqrt(2);