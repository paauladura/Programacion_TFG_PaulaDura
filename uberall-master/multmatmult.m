function A = multmatmult(varargin)

% MULTMATMULT mulitplicates multiple matrices from the left side. The
% matrices must be given in multmat-style, where each line of a matrix
% represents a rotation matrix at a certain time point.  
% Additionally, you can specify that a matrix should be transposed by
% giving a parameter 't' directly after the matrix (see example below).
% 
% IN:
%    A ........ first matrix                                        [m, 9]  
%    B ........ second matrix                                       [m, 9] 
%    't' ...... specified after a matrix, indicates that this matrix should
%               be transposed.
%    (C, D, E ....... optional, third, fourth, fifth matrix         [m, 9])
% OUT:
%    A ........ multiplication of A * B * ...                       [m, 9]  
% EXAMPLE:
%    F = multmatmult(A, 't', B, C, D, 't', E) % performs the 3x3 equivalent
%                                             % of F = A' * B * C * D' * E
% USES:
%    multmat

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias Roth (MR), GI, Uni Stuttgart               
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2016-02-01: MR, brush-up help text
%    2016-01-25: MR, initial version
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

if nargin <= 1
  error('not enough input parameters!');
end

% initializations
Atrans = false; % transpose flags
Btrans = false;
cnt_start = 2;  % start value of loop counter

% get matrix A and eventually a set transpose flag
if isnumeric(varargin{1}) % must be matrix
    A = varargin{1}; 
else
    error('first input parameter must be a matrix!');
end

if nargin >= 2
  if ischar(varargin{2});
    if lower(varargin{2}) == 't'
        Atrans = true;
        cnt_start = 3;
    else
      error('to indicate a transposed matrix, only ''T'' or ''t'' is allowed!');
    end
  end
end
  
% the loop
cnt = cnt_start;
while cnt <= nargin
  % get matrix B and eventually a set transpose flag
  if isnumeric(varargin{cnt}) % must be matrix
    B = varargin{cnt};
  else
    error('input parameter %d must be a matrix!', cnt);
  end
  
  if nargin >= cnt + 1
    if ischar(varargin{cnt + 1});
      if lower(varargin{cnt + 1}) == 't'
        Btrans = true;
        cnt = cnt + 1; % increase counter to skip transpose flag
      else
        error('to indicate a transposed matrix, only ''T'' or ''t'' is allowed!');
      end
    end
  end
  
  % do the multiplication
  A = multmat(A, B, Atrans, Btrans);
  
  % reset the transpose flags
  Atrans = false;
  Btrans = false;
  % increase counter
  cnt = cnt + 1;
end


