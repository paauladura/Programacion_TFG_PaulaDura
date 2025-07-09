function R = R2(alpha, matrix_type)

% 3D rotation around Y-axis:
%   - mathematically negative (clockwise) defined rotation if applied to
%     points (moving within the same coordinate system),
%   - mathematically positive (counterclockwise) defined rotation of points
%     if the rotation is applied to the coordinate system.
%
% IN:
%   alpha ......... rotation angle in [rad]
%   matrix_type ... optional, default behaviour:
%                     if alpha scalar: 3 x 3 rotation matrix
%                     if alpha vector: multmat-style rotation matrix
%                   by setting matrix_type = 'multmat', you can force
%                   multmat-style.
% OUT:
%   R ............. rotation matrix, dependent on your input either as
%                   3 x 3 matrix (for scalars) or in multmat-style (for
%                   vectors or by setting matrix_type to 'multmat').

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


%% check the input
if min(size(alpha)) > 1
  error('TAIHEN: x needs to be a scalar or vector!');
end

if ~exist('matrix_type', 'var')
  if max(size(alpha)) == 1
    mt = '3x3';
  else
    mt = 'multmat';
  end  
else
  switch lower(matrix_type)
    case '3x3'
      if max(size(alpha)) == 1
        mt = '3x3';
      else 
        error('TAIHEN: x must be a scalar if you request 3x3 output!');
      end 
    case 'multmat'
      mt = 'multmat';
    otherwise
      error('TAIHEN: unknown matrix_type specified!');
  end
end

%% prepare the rotation matrix
switch mt
  case 'multmat'
    alpha = standing(alpha); % x must be a standing vector!
    m = length(alpha);
    R = [cos(alpha) zeros(m, 1) -sin(alpha) zeros(m, 1) ones(m, 1) zeros(m, 1) sin(alpha) zeros(m, 1) cos(alpha)];  
  otherwise % 3x3 
    R = [cos(alpha) 0 -sin(alpha); 0 1 0; sin(alpha) 0 cos(alpha)];  
end
