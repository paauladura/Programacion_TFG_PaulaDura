function [f,th,lam] = regionalshs(klm, gridext, varargin)

% REGIONALSHS spherical harmonic synthesis on any regular regional grid
%
% HOW          f = regionalshs(klm, gridext)
%     [f,th,lam] = regionalshs(klm, gridext, 'OPTION', VALUE)
%
% IN:
%     klm ........ set of SH coefficients, either in SC/CS formats 
%     gridext .... extent of the grid in radians.
%                  [co-lat_min co-lat_max long_min long_max] where 
%                  co-lat(itude) long(itude)
% OPTIONAL:
%    'lmax' ...... maximum degree/order (default: determined from klm)
%    'quant' ..... optional string argument, defining the field quantity:
%                  - 'potential' ........ (default), potential [m^2/s^2],
%                  - 'tr' ............... grav. disturbance, 1st rad. derivative [mGal],
%                  - 'trr' .............. 2nd rad. derivative [E],
%                  - 'none' ............. coefficients define the output
%                  - 'geoid' ............ geoid height [m],
%                  - 'dg', 'gravity' .... gravity anomaly [mGal], 
%                  - 'slope' ............ size of surface gradient [arcsec], 
%                  - 'water' ............ equiv. water height [m],
%                  - 'smd' .............. surface mass density [kg/m^2]. 
%    'grid' ...... optional string argument, defining the grid:
%                  - 'pole', 'mesh' ..... (default), equi-angular (n+1)*(2n+1), 
%                                         including poles and Greenwich meridian.
%                  - 'block', 'cell' .... equi-angular block midpoints. n*2n
%    'dlat' ...... grid size in the latitude direction in radians (def.: pi/(2*lmax))
%    'dlong' ..... grid size in the longitude direction in radians (def.: pi/(2*lmax))
%                   If only one of 'dlat' or 'dlong' is given then that
%                   value is taken as the grid size.
%    'height' .... (default: 0), height above Earth mean radius [m].
%    'sub_normal'. if set, subtracts normal gravity field of the reference 
%                  ellipsoid (default: true)
%    'ellipsoid' . Reference ellipsoid
%                  - 'wgs84'
%                  - 'grs80'
%                  - 'he'
%    'legendre' .. which algorithm to use for Legendre functions:
%                  - 'plm' ... unstable (for d/o > ~1800) PLM
%                  - 'mex' ... X-number stabilized LEGENDRE_MEX
%
% OUT:
%    f ........... the global field
%    th .......... vector of co-latitudes [rad] 
%    lam ......... vector of longitudes [rad]
%----------------------------------------------------------------------------
% uses CS2SC, EIGENGRAV, PLM, NORMALKLM
% UBERALL/TWAITBAR
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU Munich
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
%----------------------------------------------------------------------------
% revision history
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users; avoidable via getopt.m)
%   2014-03-09: BD, Initial version with a lot of code inherited from GSHS
%   2015-05-26: BD, Complete rewrite of the code with the GETOPT style inputs
%                   also lots of code inherited from GSHS_
%----------------------------------------------------------------------------
% remarks
% This function is suitable only in cases, where the synthesis needs to be 
% only on a small region with regular spacing in co-latitude and longitude
% directions. 'Neumann' or 'Gauss' grids have been left as they do not make 
% sense for a regional grid.
%
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

% Primary check
if nargin < 2, error('Insufficient input arguments.'), end


%% define defaults and check optional parameters
% A lot of checking is done by EIGENGRAV as well.
defaultparams = {'lmax', inf;
                 'quant', 'potential';
                 'grid', 'mesh';
                 'dlat', inf;
                 'dlong', inf;
                 'height', 0;
                 'sub_normal', true;
                 'ellipsoid', 'wgs84';
                 'legendre', 'plm';
                 'waitbar', false}; 
params = getopt(defaultparams, false, varargin);  

% check Legendre function algorithm
plm_func = false;
switch params.legendre
    case 'mex'
        if exist('Legendre_mex', 'file') == 3; % check if compiled Legendre_mex exists
            plm_func = @Legendre_mex;
        else
            warning('Legendre_mex is not compiled! Falling back to plm...');
            plm_func = @plm;
        end
    case 'plm'
        plm_func = @plm;
    otherwise
        error('Unknown Legendre function.');
end

% field size determination, rearrange field and subtract reference field
[rw,cols] = size(klm);
if rw == cols				% klm is in CS-format
   klm = cs2sc(klm);	% convert to SC-format
elseif cols ~= 2*rw-1		% klm is in SC-format already
   error('Input "field" not in cs or sc format');
end

lmax = params.lmax;
if lmax > (rw - 1)     % desired lmax is bigger than what the field provides
    lmax = rw - 1;     % adjust lmax 
elseif lmax < (rw - 1) % if lmax is smaller than what the field provides
    klm = klm(1:(lmax+1), (rw - lmax):(rw + lmax)); % adjust field
end

% -------------------------------------------------------------------------
% Grid definition.
% -------------------------------------------------------------------------
if ~ischar(params.grid)
    error('grid argument must be string')
end 
thmin   = gridext(1); thmax   = gridext(2);
lammin  = gridext(3); lammax  = gridext(4);

% Grid size in latitude direction
if isinf(params.dlat)
    if isinf(params.dlong)
        dt = pi/lmax;
    else
        dt = params.dlong;
    end
else
    dt = params.dlat;
end

% Grid size in longitude direction
if isinf(params.dlong)
    dlam = dt;
else
    dlam = params.dlong;
end

switch lower(params.grid)
    case {'pole', 'mesh'}
        th   = (thmin:dt:thmax)';
        lam  = (lammin:dlam:lammax);
    case {'block', 'cell'}
        th   = (thmin+dt/2:dt:thmax)';
        lam  = (lammin+dlam/2:dlam:lammax);
    otherwise
        error('What type of grid do you want?')
end
nlat = length(th);
nlon = length(lam);

% -------------------------------------------------------------------------
% Preprocessing on the coefficients: 
%    - subtract reference field (if jflag is set)
%    - specific transfer
%    - upward continuation
% -------------------------------------------------------------------------
if params.sub_normal, klm = klm - cs2sc(normalklm(lmax), params.ellipsoid); 
    if params.display_warning == true
        warning('A reference field (WGS84) is removed from your coefficients')
    end
end

l       = standing(0:lmax);
transf  = eigengrav(l, params.quant, params.height);
klm     = diag(transf) * klm;

a = []; a(nlat,lmax+1) = 0;
b = []; b(nlat,lmax+1) = 0;
hwb = twaitbar('init', [], 'Regional spherical harmonic synthesis');  % initialize the waitbar

%-------------------------
% Treat m = 0 separately
%-------------------------
p = plm_func(0:lmax,0,th);
a(:,1) = p * klm(1:lmax+1,lmax+1);

%----------------------
% Do loops over orders
%----------------------
hwb = twaitbar(1 / (lmax + 1), hwb); % update the waitbar
for m = 1:lmax
    p        = plm_func(m:lmax,m,th);
    a(:,m+1) = p * klm(m+1:lmax+1,lmax+1+m);
    b(:,m+1) = p * klm(m+1:lmax+1,lmax+1-m);
    hwb      = twaitbar((m + 1) / (lmax + 1), hwb); % update the waitbar
end
twaitbar('close', hwb); % finalize the waitbar

clear klm

% -------------------------------------------------------------------------
% The second synthesis step consists of an inverse Fourier transformation
% over the rows of a and b. 
% In case of 'block', the spectrum has to be shifted first.
% -------------------------------------------------------------------------
if strcmp(params.grid, 'block') || strcmp(params.grid, 'cell') 
   m      = 0:lmax;
   cshift = ones(nlat, 1) * cos(m*pi/2/(pi/dlam));	% cshift/sshift describe the 
   sshift = ones(nlat, 1) * sin(m*pi/2/(pi/dlam));	% half-blocksize lambda shift.
   atemp  =  cshift.*a + sshift.*b;
   b      = -sshift.*a + cshift.*b;
   a      = atemp;
end

c = cos((0:lmax)'*lam);
s = sin((0:lmax)'*lam);

f = a*c + b*s;
