function c = real2cpxsh(r,lmax)

% REAL2CPXSH converts the real spherical harmonics to complex spherical
% harmonics. Use this function ony if the coefficients are normalized in
% the geodetic convention. 
%
% Klm   = (-1)^m (Clm - iSlm)/sqrt(2), m~=0
% Kl,-m = (Clm + iSlm)/sqrt(2), m~=0
% Kl0   = Cl0
%
% IN: 
%    r ....... Real spherical harmonic coefficients in CS/SC/look-up-table
%              formats
%    lmax .... Maximum degree of spherical harmonic expansion.
%
% OUT:
%    c ....... Complex spherical harmonic coefficients in CS-format
%
% USES
%    clm2sc, cs2sc, sc2cs, cssc2clm
%
% SEE ALSO:
%    cpx2realsh, real2cpxmat, real2cpxcov

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2021-04-12: MA, remove revision statement on deprecated function 
%    2010-12-04: BD, initial version
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

if nargin < 2
    error('Insufficient inputs')
end

[m,n] = size(r);

if m==(lmax^2+3*lmax+2)/2 && n==4
    [t,ind] = sort(r(:,2));
    r = r(ind,:);
    otyp = 1;
elseif m==n && m==lmax+1
    r = cssc2clm(r,lmax);
    otyp = 3;
elseif m==lmax+1 && n==2*lmax+1
    r = cssc2clm(r,lmax);
    otyp = 2;
else
    error('Input SH coefficients are not arranged in required format. Check input.')
end

c = r;
c(lmax+2:end,3) = (-1).^c(lmax+2:end,2).*(r(lmax+2:end,3) - 1i*r(lmax+2:end,4))/sqrt(2);
c(lmax+2:end,4) = (r(lmax+2:end,3) + 1i*r(lmax+2:end,4))/sqrt(2);

if otyp==2
    r = clm2sc([c(:,1:2) real(c(:,3:4))]);
    t = clm2sc([c(:,1:2) imag(c(:,3:4))]);
    
    c = r + 1i*t;
elseif otyp==3
    r = sc2cs(clm2sc([c(:,1:2) real(c(:,3:4))]));
    t = sc2cs(clm2sc([c(:,1:2) imag(c(:,3:4))]));

    c = r + 1i*t;
end
