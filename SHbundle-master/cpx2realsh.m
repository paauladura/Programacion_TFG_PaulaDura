function r = cpx2realsh(c,lmax)

% CPX2REALSH converts complex spherical harmonic coefficients to real
% spherical harmonic coefficients. Use this function only with SH
% coefficients normalized using the geodetic conventions.
%
% Klm   = (-1)^m (Clm - iSlm)/sqrt(2), m ~= 0
% Kl,-m = (Clm + iSlm)/sqrt(2), m~=0
% Kl0   = Cl0
%
% IN:
%    c ..... complex spherical harmonic coefficients in CS/SC/look-up-table
%            formats
%    lmax... Maximum degree of spherical harmonic expansion
%
% OUT: 
%    r ..... real spherical harmonic coefficients
%
% USES:
%    clm2sc, cs2sc, sc2cs, cssc2clm
%
% SEE ALSO:
%    cpx2realcov, cpx2realmat, real2cpxsh

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2021-04-12: MA, remove revision statement on deprecated function 
%    2011-05-31: BD, initial version
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

if nargin < 2, error('Insufficient input arguments'), end

[m,n] = size(c);

if m==(lmax^2+3*lmax+2)/2 && n==4
    [~, ind] = sort(c(:,2));
    c       = c(ind,:);
    otyp    = 1;
elseif m==n && m==lmax+1
    rc        = cssc2clm(real(c),lmax);
    ic        = cssc2clm(imag(c),lmax);
    c         = rc;
    c(:,3:4)  = rc(:,3:4) + 1i*ic(:,3:4);
    otyp      = 3;
elseif m==lmax+1 && n==2*lmax+1
    rc        = cssc2clm(real(c),lmax);
    ic        = cssc2clm(imag(c),lmax);
    c         = rc;
    c(:,3:4)  = rc(:,3:4) + 1i*ic(:,3:4);
    otyp = 2;
else
    error('Input SH coefficients are not arranged in required format. Check input.')
end

r = c;

r(lmax+2:end,3) = real(c(lmax+2:end,4))*sqrt(2);
r(lmax+2:end,4) = imag(c(lmax+2:end,4))*sqrt(2);

if otyp==2
    r = clm2sc(r);
    r = real(r);
elseif otyp==3
    r = sc2cs(clm2sc(r));
    r = real(r);
end
