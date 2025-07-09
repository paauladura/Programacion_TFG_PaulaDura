function [f, theRAD, lamRAD] = gshs_(field, varargin)
% GSHS_ global spherical harmonic synthesis 
%
% f = gshs_(field)
%
% IN:
%    field ....... set of SH coefficients, either in |c\s| or /s|c\ format
% OPTIONAL:
%    'max_lm' .... maximum degree/order (default: determined from field)
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
%    'sub_WGS84' . if set, subtracts reference ellipsoid WGS84 (default: true)
%    'grid' ...... optional string argument, defining the grid:
%                  - 'pole', 'mesh' ..... (default), equi-angular (n+1)*2n, including 
%                                         poles and Greenwich meridian.
%                  - 'block', 'cell' .... equi-angular block midpoints. n*2n
%                  - 'neumann', 'gauss' . Gauss-grid (n+1)*2n
%    'gridsize' .. grid size parameter n. (default: n = lmax, determined from field)
%                  #longitude samples: 2*n
%                  #latitude samples n ('blocks') or n+1.
%    'height' .... (default: 0), height above Earth mean radius [m].
%    'legendre' .. which algorithm to use for Legendre functions:
%                  - 'plm' ... (default) unstable (for d/o > ~1800) PLM
%                  - 'mex' ... X-number stabilized LEGENDRE_MEX
% OUT:
%    f ........... the global field
%    theRAD ...... vector of co-latitudes [rad] 
%    lamRAD ...... vector of longitudes [rad]
%
% EXAMPLE: see SHBUNDLE/examples/example_gshs.m
%
% USES:
%    vec2cs, cs2sc, eigengrav, plm, Legendre_mex, normalklm, 
%    UBERALL/grule, UBERALL/standing, UBERALL/isint, UBERALL/ispec,
%    UBERALL/twaitbar
% 
% SEE ALSO:
%    GSHS_GRID, GSHS_PTW, REGIONALSHS, GSHA

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias WEIGELT (MW), DoGE, UofC  
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2021-04-12: MA, remove comment/typo on deprecated function 
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users; avoidable via getopt.m)
%    2015-05-21: MR, reprogram parameter interface, code speed up, revise 
%                    help text
%    2014-03-09: BD, changed the variable 'grid' to 'grd' as 'grid' conflicts
%                    with the function 'grid'
%    2014-01-14: MR, brush up help text, beautify code, exchange deprecated
%                    'isstr' by 'ischar'
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments/removing of smoothing option
%    2013-01-23: MA, output in radian
%    2013-01-18: MA, replacement: isscal -> isscalar
%    1998-08-28: NS, brushing up and redefinition of checks
%    1994-??-??: NS, initial version
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

%% define defaults and check optional parameters
% A lot of checking is done by EIGENGRAV as well.
defaultparams = {'max_lm', inf;
                 'quant', 'potential';
                 'grid', 'mesh';
                 'gridsize', inf;
                 'height', 0;
                 'sub_wgs84', true;
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
[row, col] = size(field);

if col == row              % field in |C\S|-format 
    field = cs2sc(field);  % convert to /S|C\-format
elseif col ~= 2 * row - 1
   error('Input "field" not in cs or sc format');
% else: field is already in /S|C\-format
end

lmax = params.max_lm;
if lmax > (row - 1)     % desired max_lm is bigger than what the field provides
    lmax = row - 1;     % adjust max_lm 
elseif lmax < (row - 1) % if max_lm is smaller than what the field provides
    field = field(1:lmax+1,(row - lmax):(row + lmax)); % adjust field
% else: everything is ok
end

if isinf(params.gridsize) % no gridsize specified --> determine from field
    n = lmax;
elseif (isscalar(params.gridsize) && isint(params.gridsize)) 
    n = params.gridsize;  % otherwise take from parameter
else
    error('gridsize must be integer & scalar')
end

% -------------------------------------------------------------------------
% Grid definition.
% -------------------------------------------------------------------------
if ~ischar(params.grid)
    error('grid argument must be string')
end 
dt = pi / n;
switch lower(params.grid)
    case {'pole', 'mesh'}
        theRAD = (0:dt:pi)';
        lamRAD = (0:dt:2*pi-dt);
    case {'block', 'cell'}
        theRAD = (dt/2:dt:pi)';
        lamRAD = (dt/2:dt:2*pi);
    case {'neumann', 'gauss'}
        [gx, ~] = grule(n+1);
        theRAD = flipud(acos(standing(gx)));
        lamRAD = (0:dt:2*pi-dt);
    otherwise
        error('What type of grid do you want?')
end
nlat = length(theRAD);
nlon = length(lamRAD);

% -------------------------------------------------------------------------
% Preprocessing on the coefficients: 
%    - subtract reference field (if params.sub_wgs84 is set)
%    - specific transfer
%    - upward continuation
% -------------------------------------------------------------------------
if params.sub_wgs84, 
    field = field - cs2sc(normalklm(lmax, 'wgs84')); 
    if params.display_warning == true
        warning('A reference field (WGS84) is removed from your coefficients')
    end
end

l = standing(0:lmax);
transf = eigengrav(l, params.quant, params.height);
field  = field .* (transf * ones(1, 2 * lmax + 1));
% -------------------------------------------------------------------------
% Size declarations and start the waitbar:
% Note that the definition of dlam causes straight zero-padding in case N > L.
% When N < L, there will be zero-padding up to the smallest integer multiple
% of N larger than L. After the Fourier transformation (see below), the
% proper samples have to be picked out, with stepsize dlam.
% -------------------------------------------------------------------------
dlam   = ceil(lmax / n);				% longitude step size
abcols = dlam * n + 1;				    % # columns required in A and B
a = []; a(nlat, abcols) = 0;            % faster than a = zeros(nlat, abscols)
b = []; b(nlat, abcols) = 0;            
hwb    = twaitbar('init', [], 'gshs');  % initialize the waitbar

% treat m = 0 separately:
a(:, 1) = plm_func(0:lmax, 0, theRAD) * field(1:lmax+1, lmax+1); 
% b is already filled by '0' ... unnecessary: b(:, m+1) = zeros(nlat, 1);

hwb = twaitbar(1 / (lmax + 1), hwb); % update the waitbar
for m = 1:lmax % rest of order m
   p         = plm_func(m:lmax, m, theRAD);
   a(:, m+1) = p * field(m+1:lmax+1, lmax+1+m); % plm * c
   b(:, m+1) = p * field(m+1:lmax+1, lmax+1-m); % plm * s
   hwb       = twaitbar((m + 1) / (lmax + 1), hwb); % update the waitbar
end
twaitbar('close', hwb); % finalize the waitbar

clear field

% -------------------------------------------------------------------------
% The second synthesis step consists of an inverse Fourier transformation
% over the rows of a and b. 
% In case of 'block', the spectrum has to be shifted first.
% When no zero-padding has been applied, the last b-coefficients must be set to
% zero. Otherwise the number of longitude samples will be 2N+1 instead of 2N.
% For N=L this corresponds to setting SLL=0!
% -------------------------------------------------------------------------
if strcmp(params.grid, 'block') || strcmp(params.grid, 'cell') 
   m      = 0:abcols-1;
   cshift = ones(nlat, 1) * cos(m*pi/2/n);	% cshift/sshift describe the 
   sshift = ones(nlat, 1) * sin(m*pi/2/n);	% half-blocksize lambda shift.
   atemp  =  cshift .* a + sshift .* b;
   b      = -sshift .* a + cshift .* b;
   a      = atemp;
end

if rem(n, lmax) == 0				% Case without zero-padding.
   b(:, abcols) = zeros(nlat, 1);	
end

f = ispec(a', b')';
if dlam > 1, f = f(:, 1:dlam:dlam*nlon); end
