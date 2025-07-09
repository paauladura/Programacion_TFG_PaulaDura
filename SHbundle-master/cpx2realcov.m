function C = cpx2realcov(lmax)

% CPX2REALCOV generates a matrix the size of a spectral covariance matrix
% for a given spherical harmonic degree. This can be used for converting
% complex spherical harmonic covariances to real spherical harmonic
% covariances. This will find use mainly in rotation of spherical
% harmonics. The output is in degree-leading format.
%
% IN:
%    lmax .. Maximum degree of spherical harmonic expansion
%
% OUT: 
%    C ..... Matrix that converts the real spherical harmonics to complex
%
% USES: 
%    cpx2realmat
% 
% SEE ALSO:
%    cmp2realmat, cpx2realsh

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

lmax = floor(abs(lmax));

c = sum(1:lmax+1);
s = sum(1:lmax);

C    = struct('cc',zeros(c),'cs',zeros(c,s),'sc',zeros(s,c),'ss',zeros(s));
C.cc = mat2cell(C.cc,1:lmax+1,1:lmax+1);
C.cs = mat2cell(C.cs,1:lmax+1,1:lmax);
C.sc = mat2cell(C.sc,1:lmax,1:lmax+1);
C.ss = mat2cell(C.ss,1:lmax,1:lmax);

C.cc{1,1} = 1;

for k = 1:lmax
    tmp = cpx2realmat(k);
    
    C.cc{k+1,k+1} = tmp(k+1:end,k+1:end);
    C.ss{k,k}     = rot90(tmp(1:k,1:k),2); % = fliplr(flipud(tmp(1:k,1:k)));
    C.sc{k,k+1}   = flipud(tmp(1:k,k+1:end));
    C.cs{k+1,k}   = fliplr(tmp(k+1:end,1:k));
end
C = [cell2mat(C.cc) cell2mat(C.cs); cell2mat(C.sc) cell2mat(C.ss)];
