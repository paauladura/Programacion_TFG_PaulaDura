function [basin_mask, basins, uniq_bas, shp] = hydrosheds2pix(fnm, grid_size)

% HYDROSHEDS2PIX reads a hydroSHEDS shapefile, converting it into a 
% basin mask on an equiangular grid. The pixels contain the basin ID.
%
% [pixels, shp] = shp2pix(fnm)
%
% INPUT
% fnm       - Filename with the full path to the shapefile
% grid_size - Grid size of the basin mask in degrees
%
% OUTPUT
% basin_mask    - Equiangular grid with the grid cells conatining the 
%                 basin ID
% basins        - A cell array containing the unique basin ID, it's
%                 mask and its area
% uniq_bas      - ID of the catchments that were read from the shapefile
% shp           - A structure file containing the shapefile information
%                 as returned by the SHAPEREAD function
%
%-----------------------------------------------------------------------




% 
%  Project: GRACE Bundle
%  Copyright Balaji Devaraju (BD) 
%  devaraju at ife dot uni-hannover dot de
% 
%  License: GNU GPLv3 or later
%  You should have received a copy of the GNU General Public License
%  along with gracebundle;  If not, see <http://www.gnu.org/licenses/>.
%  
%  Authors: Balaji Devaraju
% 
%  Version control
%  Auhtor   YYYY:MM:DD  Comment
%   BD      2017:12:20  Initial version
%----------------------------------------------------------------------
% 


narginchk(1, 2)

if ~ischar(fnm)
    error('Input is not a string. Expecting a string argument.')
end

if nargin == 1
    grid_size = 0.5;
end


fprintf('\n\n---------------------------------------------------------------\n\t')
fprintf('Converting HydroSHEDS shapefile into a basin mask\n')
fprintf('---------------------------------------------------------------\n')


fprintf('Preliminaries ...\n');
shp = shaperead(fnm);

main_bas = vertcat(shp.MAIN_BAS);

% Generating grid
fprintf('Initialising X and Y grids ...\n');
[Xq, Yq] = meshgrid(-180+grid_size/2:grid_size:180,90-grid_size/2:-grid_size:-90);
[nlat, nlong] = size(Xq);
Xq = Xq(:);
Yq = Yq(:);


% In or out or on polygon
fprintf('Checking for pixels that lie within a basin ... 0000');
Rb = cell(length(main_bas), 2);
for k = 1:length(main_bas)
    fprintf('\b\b\b\b%04d', k)
    X = shp(k).X;
    Y = shp(k).Y;
    
    [in_, on_] = inpolygon(Xq, Yq, X, Y);
    
    Inside  = reshape(in_, nlat, nlong);
    On      = reshape(on_, nlat, nlong);
    
    Rb{k, 1} = main_bas(k);
    Rb{k, 2} = Inside+On;
end
fprintf('\b\b\b\bdone\n')

% Generating basin mask
fprintf('Generating basin mask ... 0000')
[uniq_bas, ~, idx]  = unique(main_bas);
basins              = cell(length(uniq_bas), 3);
basin_mask          = zeros(size(Xq));
keep_it             = true(length(uniq_bas), 1);
for k = 1:length(uniq_bas)
    fprintf('\b\b\b\b%04d', k)
    tmp         = find(main_bas==uniq_bas(k));
    basins{k,1} = uniq_bas(k);
    basins{k,2} = 0;
    basins{k,3} = 0;
    for j = 1:length(tmp)
        basins{k,2} = basins{k,2} + Rb{tmp(j), 2};
        basins{k,3} = basins{k,3} + shp(tmp(j)).SUB_AREA;
    end
    indx             = (basins{k,2} ~= 0 );
    basins{k,2}      = double(indx);
    basin_mask(indx) = uniq_bas(k);
    if ~any(indx(:))
        keep_it(k) = false;
    end
end
fprintf('\b\b\b\bdone\n')
basins      = basins(keep_it, :);
uniq_bas    = uniq_bas(keep_it);


basin_mask = reshape(basin_mask, nlat, nlong);



% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
