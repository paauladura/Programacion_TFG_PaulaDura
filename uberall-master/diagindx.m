function [ indx ] = diagindx(sz)

% DIAGINDX provides the indices of the diagonal elements of a given square matrix.
% This function is useful for replacing the values of the diagonal of a given
% matrix.
%
% IN:
%    sz ....... size of the square matrix
%
% OUT:
%    indx ..... indices of the diagonal of the sqaure matrix
%
% USES:
%    uberall/diagindx

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-03-03: BD, initial version
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

r 	= (0:sz-1); % rows
c 	= (1:sz);   % columns
indx 	= r*sz + c; % indices
indx 	= indx(:);
