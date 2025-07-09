function [xy,theRAD,lamRAD] = cindxcoord(cindx,grd,ldir,indx)

% CINDXCOORD finds the co-ordinates of various catchments and provides
% them in a cell array, with each cell representing one catchmenta. The 
% output is in radians.
%
% xy = cindxcoord(cindx,grd,ldir)
% xy = cindxcoord(cindx,grd,ldir,indx)
%
% INPUT
% cindx     -   A matrix providing the indices of the catchments.
% grd       -   Specify the type of the grid: 'block'/'cell', 'mesh'/'pole', 
%               'gauss'/'neumann'.
% ldir      -   '1' indicates that Longitude is arranged 0-360 and '0' 
%               indicates -180 to 180 arrangement.
% indx      -   Indices of the specific catchments whose coordinates are 
%               required.
% OUTPUT
% xy        -   Cell array containing cells for each catchment and their
%               corresponding co-ordinates in the grid. [radians]
%               {index [x y]}
% theRAD    -   Co-latitude of the grid.
% lamRAD    -   Longitude of the grid.
%--------------------------------------------------------------------------

% Created on: 25 March 2008, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

% Checking inputs
if nargin < 3
    error('Insufficient input arguments')
elseif nargin == 3
    indx = unique(cindx(:));
end


% ------------------------
% Grid definition.
% ------------------------
n = size(cindx,1);
if (mod(n,2)~=0), n = n - 1; end

dt = pi/n;
if strcmp(grd,'pole') | strcmp(grd,'mesh')
   theRAD   = (0:dt:pi)';
   lamRAD   = (0:dt:2*pi-dt);
elseif strcmp(grd,'block') | strcmp(grd,'cell')
   theRAD   = (dt/2:dt:pi)';
   lamRAD   = (dt/2:dt:2*pi);
elseif strcmp(grd,'neumann') | strcmp(grd,'gauss')
   [gx,gw]  = grule(n+1);
   theRAD   = flipud(acos(standing(gx)));
   lamRAD   = (0:dt:2*pi-dt);
else
   error('What type of grid do you want?')
end

ldir = logical(ldir);
if lam
    lamRAD(lamRAD > pi) = lamRAD(lamRAD > pi) - 2*pi;
    lamRAD              = sort(lamRAD);
end

nth     = length(theRAD);
nlam    = length(lamRAD);

X = ones(nth,1) * lamRAD;
Y = theRAD * ones(1,nlam);

X = X(:);
Y = Y(:);
c = cindx(:);

xy = cell(size(indx(:),1),2); % Initializing cell array
 
for i = 1:size(indx(:),1)
    xy{i,1} = indx(i); % Writing index
    cnt = (c==indx(i));
    xy{i,2} = [X(cnt) Y(cnt)]; % Writing co-ordinates
end
