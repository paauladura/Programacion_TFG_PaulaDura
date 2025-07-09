function f = ispec2(cc,cs,sc,ss)

% ISPEC2 returns 2D function from real Fourier coefficients
%     
% IN:
%    cc ...... cos-cos coefficient
%    cs ...... cos-sin coefficient
%    sc ...... sin-cos coefficient
%    ss ...... sin-sin coefficient
%
% OUT:
%    f ....... 2D function (matrix)
% 
% USES: 
%    ispec, spec
%
% SEE ALSO:
%    SPEC2, FFT2, SPEC, ISPEC
%
% REMARKS:
%    Defining whether size(f,i) will be even or odd is based on condition
%    (sine-coefficient < 1e-10) in ISPEC. This is aliasing at Nyquist freq.
%    This might produce wrong results after low-pass filtering.

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
c = ispec(cc',sc')';
s = ispec(cs',ss')';
f = ispec(c,s);

