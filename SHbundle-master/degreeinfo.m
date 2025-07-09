function d = degreeinfo(field,cum,quant,h)

% DEGREEINFO transforms spherical harmonic information by degree. Input can be 
% SH-coefficients, their SD's, signal RMS values or error RMS's.
% The output is a vector of "degree information", corresponding to l=2:L,
% valid for some gravitational functional.  
% Output may be chosen cumulative. If not, output is in RMS-sense.
%
% IN: 
%    field .. either a matrix of (SD's of) spherical harmonic coefficients 
%             (SC-triangle or CS-square format), or a vector of degree RMS's,
%             in which case it corresponds to l=0:L
%    cum .... accumulate the degree info, eg. commission error, or not.
%             cum = 0 for RMS values (default)
%             cum = 1 for cumulatives
%             cum = 2 for root median square values with division by 2*l+1
%             cum = 3 for root median square without division
%             cum = 4 for median of absolute only
%             cum = 5 for variances
%             cum = 6 for RMS ignoring 0 elements
%    quant .. optional string argument, defining the field quantity:
%             'none' (default), 'geoid', 'dg' or 'gravity' (grav. anomaly),
%             'tr' (grav. disturbance), 'trr' (d^2/dr^2), or 'potential'
%    h ...... height [m]
%             default: 0
% OUT:
%    d ...... vector of SH-spectral "information" by degree
%
% USES:
%    cs2sc, eigengrav, sc2cs, uberall/standing

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%   2013-02-13: MR, change function names, brush up comments
%   2013-01-30: MA, comments/removing of smoothing option
%   1995-05-09: NS, initial version
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

% Defaults and further checks. A lot of checking will be done by EIGENGRAV.
if nargin < 4, h     = 0;      end
if nargin < 3, quant = 'none'; end
if nargin < 2, cum   = 0;      end
if nargin < 1, error('Input please.'), end
if isempty(h),     h     = 0;      end
if isempty(quant), quant = 'none'; end
if isempty(cum),   cum   = 0;      end

% Checks and degree-variance determination:
[rows,cols] = size(field);
if isvector(field)				    % degrms's as input
   field = standing(field);
   lmax  = length(field)-1;
   l     = standing(0:lmax);
   dv    = (2*l+1) .* (field.^2);	% degrms turned into degvar
   dv(1:2) = [];					% eliminate l = 0 and 1
   
   % Computation of transfer (lambda_l) vector:
   l(1:2) = [];					    % eliminate l = 0 and 1
   tf = eigengrav(l,quant,h);
   
   % Actual computation:
   if cum == 2
       error('Median for input degrms not possible.')
   elseif cum == 1
       d = sqrt(cumsum(dv .* tf.^2));
   else
       d = sqrt(dv ./ (2*l+1) .* tf.^2);
   end
   return
elseif rows == cols				    % field is in CS-format
    lmax  = rows - 1;
    field = cs2sc(field,0);	        % degree variances
    l     = standing(0:lmax);
elseif cols-2*rows == -1			% field is in SC-format already
    lmax  = rows - 1;
    field = cs2sc(sc2cs(field),0);
    l     = standing(0:lmax);
else
    error('Check format of FIELD.')
end
dv = field;

% Computation of sum or median
if cum == 4
    helper = zeros(size(dv,1),1);
    for ridx = 0:lmax
        helper(ridx+1) = median(dv(ridx+1,lmax+1-ridx:lmax+1+ridx));
    end
    dv = helper;
elseif cum == 2 || cum == 3
    helper = zeros(size(dv,1),1);
    for ridx = 0:lmax
        helper(ridx+1) = median(dv(ridx+1,lmax+1-ridx:lmax+1+ridx).^2);
    end
    dv = helper;
else
    dv = sum(dv.^2,2);
end

% Computation of transfer (lambda_l) vector:
dv(1:2) = [];					% eliminate l = 0 and 1
l(1:2)  = [];					% eliminate l = 0 and 1
tf = eigengrav(l,quant,h);

% Actual computation:
if cum == 6
    field = cs2sc(sc2cs(field),NaN);
    d = sqrt(dv ./ (2*(l-sum((field(3:end,:)==0),2))+1) .* tf.^2); 
elseif cum == 5
    d = sqrt(dv.*tf.^2);
elseif cum == 4
    d = abs(dv).*tf;
elseif cum == 3
    d = sqrt(dv.*tf.^2);
elseif cum == 2
    d = sqrt(dv ./ (2*l+1) .* tf.^2);
elseif cum == 1
    d = sqrt(cumsum(dv .* tf.^2));
else
    d = sqrt(dv ./ (2*l+1) .* tf.^2);
end

