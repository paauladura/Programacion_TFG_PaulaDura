function [a,b] = spec(y)

% SPEC(Y) returns the cosine and sine spectra of Y.
%
% OUT:
%     a ...... cosine coefficients
%     b ...... sine coefficients
%              The a and b coefficients are defined by:
%              y(t) = a_0 + SUM_(i=1)^n2 a_i*cos(iwt) + b_i*sin(iwt)
%   
%              with w = ground-frequency and n2 half the number of samples (+1).
%              Note that no factor 2 appears in front of the sum.
% 
%             [A,B] = SPEC(Y) will return the cosine and sine spectra separately.
%              S = SPEC(Y) or just SPEC(Y) will return a matrix [A B].
%              If Y is a matrix, spectra are computed columnwise.
% 
% SEE ALSO:
%    ISPEC, FFT, etc.

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1994-04-27: NS, initial version
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

if min(size(y)) == 1
   y = y(:); 				                % if vector, put y upright.
end
[n m] = size(y);

sy = fft(y);
n2 = fix((n+2)/2);
sy = sy(1:n2,:);

a = (sy + conj(sy)) / n;
a(1,:) = a(1,:) / 2;
if rem(n,2) == 0, a(n2,:) = a(n2,:)/2; end	% eliminate the 100% aliasing
b = 1i * (sy - conj(sy)) / n;

if nargout < 2
   a = [a b];
end
