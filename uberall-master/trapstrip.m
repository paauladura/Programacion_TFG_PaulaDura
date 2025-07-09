function b = trapstrip(a)

% TRAPSTRIP(A) strips the upper and lower left corners from matrix A
% (standing) or the left and right upper corners (A lying). 
% The corner elements are replaced with zeros. A trapezium shape remains.
%
% SEE ALSO:
%    TRIL, TRIU

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1994-03-08: NS, initial version
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

[rows,cols] = size(a);

if cols > rows
   b = rot90(tril(flipud(tril(fliplr(a')))),2)';
else
   b = rot90(tril(flipud(tril(fliplr(a)))),2);
end;
