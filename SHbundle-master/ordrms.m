function orms = ordrms(field)

% ORDRMS(FIELD) computes the order RMS of a set of spherical
% harmonic coefficients. It represents the average size
% of a single spherical harmonic coefficient (or its standard deviation).
%
% IN:
%    field ..... must be in SC-triangular or in CS-square format
%
% OUT:
%    orms ...... order RMS of a set of spherical harmonic coefficients
%
% USES:
%    cs2sc

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1994-12-07: NS, the number of terms per order m not defined anymore as L-m+1,
%                    but as the number of elements larger than the background.
%    1994-10-26: NS, initial version
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

[rows, cols] = size(field);
% lmax        = rows -1; 
if rows == cols
   field = cs2sc(field);
elseif cols ~= 2*rows - 1
   error('Input FIELD not in required format!')
end

% order = -lmax:lmax;
% terms = lmax + 1 - abs(order);
backv = field(1);			% background value
terms = sum(abs(field) > backv);	% #terms per m larger than background
orms  = sqrt(sum(field.^2)./terms);
orms(~terms) = NaN * ones(sum(~terms), 1); % replace Inf's by NaN's (for plotting)
