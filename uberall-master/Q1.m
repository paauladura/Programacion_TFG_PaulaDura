function Q = Q1(matrix_type, len)

% 3D mirroring at X-axis.
%
% IN:
%   matrix_type ... optional,
%                   '3x3' ....... 3 x 3 mirroring matrix,
%                   'multmat' ... multmat-style mirroring matrix, the
%                                 parameter len becomes mandatory.
%   len ........... only if matrix_type is set to 'multmat': length of
%                   vector-array.
% OUT:
%   Q ............. matrix for mirroring at X-axis, dependent on your input
%                   either as 3 x 3 matrix (default) or in multmat-style
%                   (by setting matrix_type to 'multmat' and len to the
%                   length of the vector-array). 

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias Roth (MR), GI, Uni Stuttgart               
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2016-02-01: MR, brush-up help text
%    2016-01-27: MR, initial version
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


if ~exist('matrix_type', 'var') && ~exist('len', 'var')
  len = 1;
elseif ~exist('len', 'var');
  switch lower(matrix_type)
    case '3x3'
      len = 1;
    otherwise
      error('TAIHEN: specify "len" or I''m unable to generate the output!');
  end
end

%% prepare the rotation matrix
if len > 1
    Q = [-ones(len, 1) zeros(len, 3) ones(len, 1) zeros(len, 3) ones(len, 1)]; 
else % 3x3 
    Q = [-1 0 0; 0 1 0; 0 0 1];
end