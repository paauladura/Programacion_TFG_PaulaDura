function [xx,yy,zz] = xyz2mat(xyz,backval,zcol)

% XYZ2MAT(xyz,backval,zcol) converts column data XYZ into a matrix ZZ. 
%
%       It is assumed that the 1st and 2nd column contain x and y information,
%       not necessarily organized in a specific way. Gaps will be filled with 
%       the background value BACKVAL, which has default zero. If z-data is 
%       stored in a column different from the 3rd one, a specific ZCOL value 
%       must be supplied as well.
%
% REMARKS:
%     - The x and y information must be INTEGER. All the x- and y-value
%       will be shifted, such that their minimum has the value one. 
%       Corresponding matrices XX and YY can optionally be created though.
%       In that case, use [xx,yy,zz] = XYZ2MAT(xyz,backval,zcol)
%     - Double entries will be overwritten.
%     - Check "full(spconvert(xyz))" as well.
%
% USES:
%    isint

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1995-12-12: NS, No separate vectors x,y,z taken out of xyz anymore
%    1995-11-17: NS, Problem of x,y vs. i,j taken care of
%    1995-03-23: NS, vectorized mapping, by defining a pointer
%    1994-03-14: NS, initial version
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

% Defaults and checks
if nargin == 2
   zcol = 3;
elseif nargin == 1
   backval = 0;
   zcol = 3;
elseif nargin ~= 3
   error('Check the (number of) input arguments')
end
if ~all(all(isint(xyz(:,1:2)))), error('INTEGER xy-columns, please'), end


% Preliminaries
xymin = min(xyz(:,1:2));
xymax = max(xyz(:,1:2));
xmin  = xymin(1);
ymin  = xymin(2);
xmax  = xymax(1);
ymax  = xymax(2);
n     = xmax-xmin+1;
m     = ymax-ymin+1;


% Do the mapping
zz      = backval * ones(m,n);				% init
ptr     = (xyz(:,1)-xmin)*m + xyz(:,2)-ymin+1;		% pointer
zz(ptr) = xyz(:,zcol);					% very compact

% Output
if nargout == 1
   xx = zz;
elseif nargout == 3
   [xx,yy] = meshgrid(xmin:xmax,ymin:ymax);
end
