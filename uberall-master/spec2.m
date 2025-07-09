function [cc, cs, sc, ss] = spec2(f)

% SPEC2 returns 2D real Fourier coefficients
%
% HOW:
%    [cc, cs, sc, ss] = spec2(f)
%                   p = spec2(f)
%     
% IN:
%    f ....... 2D function (matrix)
%
% OUT:
%    cc ...... cos-cos coefficient
%    cs ...... cos-sin coefficient
%    sc ...... sin-cos coefficient
%    ss ...... sin-sin coefficient
%    p ....... power spectrum (sqrt)
%
% USES:
%    spec
%
% SEE ALSO:
%    ISPEC2, FFT2, SPEC, ISPEC
%
% REMARKS:
%    Fastest when size(f) is power of 2.

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1998-02-19: NS, initial version
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

% here we go
[c,s]   = spec(f);
[cc,sc] = spec(c');
[cs,ss] = spec(s');
cc=cc'; cs=cs'; sc=sc'; ss=ss';


if nargout < 2
   cc = sqrt(cc.^2+cs.^2+sc.^2+ss.^2);
end
