function out = replace(in,where,by)

% REPLACE puts a replacement whenever an element of the input vector
% fulfills some condition (denoted by index or mask)
%      
% IN:
%    in ........ input		      	- vector (or matrix)
%    where ..... index vector or logical
%    by ........ replacement 		- scalar, vector or []
%
% OUT:
%    out ....... output		     	- vector (or matrix)
%
% REMARKS:
%    - Works only on matrix if BY is scalar or empty.
%    - If BY is empty matrix, replace becomes a delete operation. In case IN is
%      a matrix, it is checked whether WHERE indicates rows or columns alone. 
%      In that case OUT will be matrix, otherwise vector. 
%
%    - perhaps not elegant, but fully vectorized
%    - not tested yet on strings
%    - nice example of index juggling
%    - could be extended for the case with matrix IN, vector BY and
%      row/column case of where

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-04-08: BD, changed short-circuit '&' to '&&'
%    2013-02-13: MR, change function names, brush up comments
%    1998-12-23: NS, V5-update: - scalar assignment implemented
%                               - check whether WHERE is a logical 
%    1997-08-08: NS, if isempty(BY), it is checked whether the delete operation
%                    can be performed on rows or columns alone -> matrix output
%    1997-06-16: NS, bug fix - in case where=1, it is identified as index now
%    1997-06-09: NS, empty matrix BY possible now
%    1996-12-06: NS, initial version
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

% diagnostics and preliminaries
narginchk(3,3);
% If WHERE is an index, transform into mask. 
if ~islogical(where)
   idx   = where;
   where = logical(zeros(size(in)));
   where(idx) = 1;
end
if ~all(size(where)==size(in)), error('IN and WHERE don''t match'), end


% here we go
nby   = length(by);
[r,c] = size(in);

if isempty(by) 				% BY is empty. Also for matrix IN

   out = in;
   cmsk = sum(where);
   rmsk = sum(where');
   if all(diff(cmsk) == 0) && all(rmsk == 0 | rmsk == c)
						% row case. OUT will be matrix
      out(sign(rmsk),:) = [];
   elseif all(diff(rmsk) == 0) && all(cmsk == 0 | cmsk == r)
						% column case. OUT will be matrix
      out(:,sign(cmsk)) = [];
   else					% general case. OUT will be vector
      out(where) = [];
   end

elseif nby == 1				% BY is scalar. Also for matrix IN

   out = in;
   out(where) = by;
else						% otherwise

   if min(size(in)) > 1
      error('non-scalar REPLACE only operates on vectors') 
   end
   isrowvec = (diff(size(in)) > 0);
   in   = in(:);
   by   = by(:);
   where = where(:);
   nin  = length(in);
   nocc = sum(where);			% #occurrences
 % first, shift the elements of IN
   out  = zeros(nin+nocc*(nby-1),1);	% initialize output vector 
   idx  = cumsum(where*(nby-1)) + (1:length(in))'; 
   out(idx) = in;
 % next, replace the indicated elements
   by   = kron(ones(nocc,1),by);		% replacement vector [by;by;by;...]
   idx  = find(where) + (nby-1)*(0:nocc-1)';
   idx  = cumsum([idx(:)'; ones(nby-1,nocc)]);
   out(idx(:)) = by;			% do the replacement
   if isrowvec, out = out'; end		% row vector if necessary

end

