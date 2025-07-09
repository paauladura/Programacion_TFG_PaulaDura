function out = shprepare(in,how,lmax)

% SHPREPARE prepares a spherical harmonic spectrum for imaging purposes.
% It is turned into SC-format, if not already.
% Returned is the INput field itself, its gain, or significant digits.
% OUTput is logarithmic.
% The background value is one (order of magnitude) less than foreground.
%
% IN:
%    in ..... SH-spectrum in CS or SC format
%    how .... flag: return field itself (how=0), its gain (how=1), or
%             significant digits (how=2)
%    lmax ... maximum degree. not necessarily the same as from input field.
%
% OUT:
%    out .... logarithmic SH-spectrum in SC-format
% 
% USES: 
%    cs2sc, sc2cs, kaula, jgm1s.mat

% REMARKS:
%    - Gain is defined w.r.t. JGM1s errors and Kaula.
%    - Significant digits are w.r.t. Kaula. (floating number)
%    - how=1 and how=2 are only relevant for SH error spectra
%    - if LMAX is chosen larger than the SH spectrum, Kaula's rule is padded
%      at the end. This implies zero gain and zero significant digit.

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1998-12-21: NS, Matlab V5 update: logical mask
%    1996-09-24: NS, - HOW-flag added: gain and significan digits option
%                    - flexible LMAX-option added
%    1996-09-24: NS, initial version
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

% init
[rows,cols] = size(in);
lmaxin      = rows-1;
if nargin < 3, lmax = lmaxin; end
if nargin < 2, how = 0; end		% return field itself (default)
if ~any(0:2==how), error('Set HOW-flag properly.'), end

% Put the INput in SC-format.
if cols == 2*rows - 1			% in SC-format 
    in = sc2cs(in);
elseif rows ~= cols
   error('Input FIELD not in required format!')
end

lvec   = 1:lmax+1;
lvecin = 1:lmaxin+1;
if lmax < lmaxin
   in = in(lvec,lvec);			% shrink IN to new LMAX
end
k   = [1e-20 1e-20 kaula(2:lmax)]';
kau = sc2cs(k*ones(1,2*lmax+1));

% Take logarithm of KAU and IN already. Subsequent calculations are additive.
in  = log10(abs(in));
kau = log10(kau);

% Create mask
mask = ones(lmax+1);                 % create mask for valid coeffs (CS)
mask(1:2,1:2) = zeros(2);            % l=0 and l=1 terms 
mask(3,2) = 0; mask(1,3) = 0;        % C21 and S21
mask = logical(mask);			% necessary for V5

% Now get OUT for the case HOW=0. Also basis for HOW = 1 or 2.
if lmax > lmaxin
   out = kau;
   out(lvecin,lvecin) = in;
else
   out = in;
end

% Further calculations for HOW = 1 or 2.
if how == 1
   load jgm1s
   jgm1s_sd = log10(jgm1s_sd);
   bv = -1;
   if lmax > 60
      ref = kau;
      ref(1:61,1:61) = jgm1s_sd;
   else
      ref = jgm1s_sd(lvec,lvec);
   end
   out = max(ref - out,0);
elseif how == 2
   bv  = -1;
   out = max(kau - out,0);
else 
   bv = min(out(mask))-1;
end

% Return OUT in SC-format with proper BV background values.
mask = ~mask;			% invert mask 
out(mask) = bv*mask(mask);
out = cs2sc(out,bv);


