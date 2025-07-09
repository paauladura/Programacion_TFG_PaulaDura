function C = cpx2realmat(l)

% CPX2REALMAT generates a matrix operator that converts complex spherical
% harmonics to real spherical harmonics.
%
% IN:
%    l .... spherical harmonic degree
%
% OUT: 
%    C .... matrix operator
%
% SEE ALSO:
%    real2cpxmat

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2010-11-29: BD, initial version
% ----------------------------------------------------------------------------
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
% ----------------------------------------------------------------------------

C = ones(2*l+1,1);
m = (-l:l)';

C(1:l) 		= 1i/sqrt(2); % Treating negative orders
C(l+1) 		= 0.5; % Treating order zero separately
C(l+2:end) 	= (-1).^(m(l+2:end))/sqrt(2);
C 			= diag(C) + fliplr(diag([-C(1:l); C(l+1:end)].*(-1).^(m)));
