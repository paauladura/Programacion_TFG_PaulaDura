function fbar = globalmean(f, typ, varargin)

% GLOBALMEAN computes the mean of a field given on a unit sphere. The
% field can either be in spherical harmonic coefficients or as a spatial
% map.
%
% fbar = globalmean(f, 'spatial', GRID)
% fbar = globalmean(f, 'spectral', OPTION, VALUE)
%
% INPUT
% f     - Field given on a sphere as a gridded map or as a spherical
%         harmonic spectrum
% typ   - Type of representation 
%           'spatial'
%           'spectral'
%
% SPATIAL
%--------
% grid  - If the representation is spatial then provide the type of grid
%         of the map. Valid values are
%                   'block'/'cell'      - equiangular block midpoints
%                   'mesh'/'pole'       - equiangular grid points
%                   'gauss'/'neumann'   - Gauss-neumann grid
%
% SPECTRAL OPTIONS
%-----------------
% 'quant' - If the representation is spectral then a gravity
%           functional can be sought. Valid values are
%               'none' [default]
%               'potential'
%               'tr'
%               'trr'
%               'none'
%               'geoid'
%               'gravity'
%               'slope'
%               'water'
%               'smd'
% 'height'- Height above earth mean radius [m] [default - 0]
% 'GM'    - Gravitational constant [default - IERS2010]
% 'ae'    - Semi-major axis of the Earth ellipsoid [default - IERS2010]
%
% OUTPUT
% fbar  - Global mean of the given field on the sphere
%-----------------------------------------------------------------------


% 
%  Project: GRACE Bundle
%  Copyright 2007-17 Balaji Devaraju (BD) 
%  devaraju at ife dot uni-hannover dot de
% 
%  License: GNU GPLv3 or later
%  You should have received a copy of the GNU General Public License
%  along with EGRAFS;  If not, see <http://www.gnu.org/licenses/>.
%  
%  Authors: Balaji Devaraju
% 
%  Version control
%  Auhtor   YYYY:MM:DD  Comment
%   BD      2017:05:16  Initial version
%----------------------------------------------------------------------


% Checking the number of input arguments
narginchk(2, Inf)

switch lower(typ)
    case { 'spatial' }
        if nargin < 3
            error('Cannot proceed without knowing the type of grid given')
        else
            n = size(f,1);
            grd = lower(varargin{1});

            if ~ischar(grd)
                error('Grid type is not a string')

            elseif strcmp(grd,'pole') || strcmp(grd,'mesh')
                if size(f, 2) ~= (2*(n-1))
                    error('Grid type of the field is not a "mesh"')
                end
                dt      = 180/(n - 1);
                theta   = (0:dt:180)';
                wt      = sin(theta) * dt;

            elseif strcmp(grd,'block') || strcmp(grd,'cell')
                if size(f, 2) ~= (2*n)
                    error('Grid type of the field is not "block/cell"')
                end
                dt      = 180/n;
                theta   = (dt/2:dt:180)';
                wt      = sin(theta) * dt;

            elseif strcmp(grd,'neumann') || strcmp(grd,'gauss')
                if size(f, 2) ~= (2*(n-1))
                    error('Grid type of the field is not "Gauss-Neumann"')
                end
                dt      = 180/(n - 1);
                [gx,wt] = grule(n);
                theta   = flipud(acos(standing(gx)))*180/pi;

            else
               error('Unknown grid type')
            end

            theta   = theta * pi/180;
            dt      = dt * pi/180;
            [X, Y]  = meshgrid(dt*ones(size(f,2), 1), wt);

            fbar    = bsxfun(@times, f, X);
            fbar    = bsxfun(@times, fbar, Y);
            fbar    = sum(fbar(:))/(4*pi);
        end
    case { 'spectral' }
        defpar = {  'quant', 'none'; ...
                    'height', 0; ...
                    'GM', 3.986004418e14; ...
                    'ae', 6378136.6};
        params = getopt(defpar, false, varargin);
        lam_ = eigengrav(0, params.quant, params.height, [params.gm, params.ae]);

        [chk, lmax] = checkshformat(f);

        switch lower(chk)
            case { 'cs' }
                fbar = lam_ * f(1);
            case { 'sc' }
                fbar = lam_ * f(1, lmax+1);
            case { 'clm_deg', 'clm_ord', 'klm' }
                fbar = lam_ * f(1, 3);
            otherwise
                error('Unknown format for storing the spherical harmonic coefficients')
            end
    otherwise
        error('Unknown representation of the field given on a sphere')
    end
            

            

% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
