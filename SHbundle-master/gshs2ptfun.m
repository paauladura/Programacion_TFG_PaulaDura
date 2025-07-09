function [b,th,lam] = gshs2ptfun(B,lmax,varargin)

% GSHS2PTFUN performs the spherical harmonic synthesis of the spectrum of a 
% two-point function, for example, covariance functions, smoothing operators 
% and such.
%
% b          = gshs2ptfun(B,lmax)
% [b,th,lam] = gshs2ptfun(B, lmax, 'location', VALUE, ...
%                                  'quant', VALUE, ...
%                                  'grid', VALUE, ...
%                                  'gridsize', VALUE, ...
%                                  'gridext', VALUE, ...
%                                  'gridcoord', VALUE, ...
%                                  'filter', VALUE, ...
%                                  'height', VALUE, ...
%                                  'ivartype', VALUE, ...
%                                  'const', VALUE)
%
% IN:
%    B ........ A [(L+1)^2 * (L+1)^2] variance-covariance matrix arranged
%               in Colombo ordering, or, a structure with two variables:
%               1. a character array with the covariance matrix file;
%               2. path where the file is stored (this order must be
% 				maintained at all times);
%               or, a character array of the filename and the path
%               together,
%               or, if only the variances are available then provide them
%               in either SC, CS, or Colombo ordering formats.
%               Use the first option only when the maximum degree of the
%               spherical harmonic expansion is [lmax <= 70], and also
%               provide the complete covariance matrix.
%    lmax ..... Maximum degree of the spherical harmonic expansion.
%    location . Co-ordinates for the point at which the two-point function
%               is sought. It must be of the form [90-lat,lam] (radians)
%                                                          -def: [pi,pi]/4
%    quant .... optional string argument, defining the field quantity:
%               'none' (coefficients define the output), 'geoid',
%               'potential', 'dg' or 'gravity' (grav. anomaly), 'tr' (grav.
%               disturbance), or 'trr' (2nd rad. derivative), 'water'
%               (water equivalent heights), 'smd' (surface mass density).
%               If in the case of cross-covariance computation quant is a
%               cell array with two-variables. One for each field.
%                                                           -def: 'none'
%    grid ..... optional string argument, defining the grid:
%               1. 'pole' or 'mesh': equi-angular (N+1)*2N, including poles
%               and Greenwich meridian.                     -def: 'pole'
%               2. 'block' or 'cell': equi-angular block mid-points.
%               3. 'neumann' or 'gauss': Gauss-grid (N+1)*2N. If you use
%               'gauss' the settings for 'length' settings in 'n' are
%               overridden and only a global map is provided.
%               4. 'points': If the value of the covariance function is 
%               required only a few select then this option can be used.
% 				Longitude is always counted from [0,2*pi]
%    gridcoord  Coordinates of the two-point function at 'loc' w.r.t. the
%               'loc'. Optional argument
%               1. 'polar'  - Polar coordinates with 'loc' as the pole of
%                             the two-point function at 'loc'.
%               2. 'geo'    - Geographical coordinates
%    gridsize . Grid size parameter in radians [scalar]. If the chosen grid
%               type is 'points' this value is ignored.
%    gridext .. Grid extent on either side of the location of the two-point 
%               function as a scalar or row-/column-vector. For example,
%               [pi/9,pi/6,pi/4,pi/3] tells us that the two-point function 
%               must be calculated upto 20 degrees to the north of the 
%               location, 30 degrees to the east, 45 degrees to the south, 
%               and 60 degrees to the west of the location. The vector
%               reads as [N E S W]. If an empty matrix is provided then a
%               default value of pi/9 (20 degrees) is used in all directions.
%               If only one value is given then the same length is used for
%               all the four directions. The parameter also accepts a string
%               value, which should be 'global'. This values indicates that 
%               the two-point function must be computed for the whole globe.
%               If the chosen grid type is 'points' then a two-column vector
%               must be provided for the points where the values are sought.
%               -----------------------------------------------------------
%               NOTE: The grid extent must be specified in the coordinate
%               chosen in 'gridcoord'
%               -----------------------------------------------------------
%    filter ... Spectrum of filter that needs to be applied on the two-point
%               function. Filter co-efficients should either be 
%               given as a [l Wl] vector (isotropic case) or as a matrix 
%               in CS- / SC- / [l m Clm Slm] in colombo ordering formats
%               (anisotropic case).
%                                                           - def: 'none'
%    height ... height for upward continuation [m]. It is also possible to
%               provide to different heights if the covariance matrix is
%               between two fields that are representative of two different
%               heights.
%                                                           - def: 0.
%    ivartype . type of variance-covariance propagation
%               1. 'full' - full spectrum of the two-point function is
%                           propagated.
%               2. 'block'- block-diagonal spectra of the two-point function
%                           is propagated.
%               3. 'diag' - only the diagonal spectrum of the two-point
%                           function is propagated.
% 	         	4. 'cs'	  - if input is in CS/SC/[l m Clm Slm] formats.	
%                                                           - def: 'diag'
%    const .... Earth constants for computing the certain isotropic transfer
%               functions.
%               1. If left empty constants of the GRS80 ellipsoid will be used
%               2. A 2 element array with GM and semi-major axis of the 
%                  ellipsoid can be given [GM ae]
%               3. A file with the pathname as a string can be given where 'GM'
%                  and 'ae' must have the same variable names.
%
%
% OUT:
%    b ....... covariance function for the region requested
%    th ...... vector of co-latitudes. If 'polar' was chosen for 'gridcoord'
%              then a vector of spherical distances from the center of the
%              two-point function is provided [rad]
%    lam ..... vector of longitudes. In case of 'polar' 'gridcoord' vector of
%              azimuths are provided [rad]
%
% EXAMPLE:
%    B  	= rand(46);
%    lmax   = 45;
%    loc	= [pi/2 pi/5];
%    grd 	= 'pole';
%    sz     = pi/360;
%    len    = [pi/9 pi/8 pi/6 pi/10];
%    fil    = sc2cs(clm2sc(cssc2clm([(0:lmax)' gaussfltr(pi/36,lmax)])));
%    h 		= 0;
%
%    [f,th,lam] = gshs2ptfunc(B, lmax, 'location', loc, 'grid', grd, ...
%                                      'gridsize', sz, 'gridext', len ...
%                                      'filter', fil, 'height', h, ...
%                                      'gridcoord', 'geo', 'ivartype', 'cs')
%
% USES:
%    GSHS_, GSHS_PTW, REGIONALSHS
%    CS2SC, eigengrav, PLM, COVORD, clm2sc, sc2cs, CSSC2CLM, 
%    FilterBundle\, uberall\standing, uberall\GRULE, 
%    
%
% SEE ALSO:
%    gshs_, gsha, gshscov
%

% Internal remarks:
% This is a complete rewrite of the GSHSCOVFN and it uses an alternative 
% algorithm, which is much faster for the case of 'full'ly populated B 
% matrix

% Algorithm: b = y*B*Y'
% Since this function computes the spatial form of the two-point function 
% only at specified locations, we first compute the spherical harmonics at
% those locations. This is usually a [(lmax+1)^2 n] with n << (lmax+1)^2.
% This we left multiply with the spectrum of the two-point function, which
% is again a (n x (lmax+1)^2) array. The last step is the spherical harmonic
% synthesis of the coefficients from y*B
% Step 1: Compute spherical harmonics for the specified locations 'y'
% Step 2: Compute the product 'y*B'
% Step 3: Perform spherical harmonic synthesis over each of the points yB*Y'

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2021-04-12: MA, remove comment/inactive code of deprecated function 
%    2015-03-29: BD, initial version
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

tic
fprintf('------------------------------------------------------------------\n')
fprintf('\t Computing spatial form of two-point functions \n')
fprintf('------------------------------------------------------------------\n')


fprintf('Checking input arguments ... ')
%-----------------------------------------------------
% Checking the input and augmenting initial values
%-----------------------------------------------------
params = struct('location', [pi/4,pi/4], ...
                'quant', 'none', ...
                'grid', 'mesh', ...
                'gridsize', pi/180, ...
                'gridext', pi/9, ...
                'gridcoord', 'geo', ...
                'filter', 'none', ...
                'height', 0, ...
                'ivartype', 'diag', ...
                'const', []);
pnames = fieldnames(params);

%% check for specified parameters    
if rem(length(varargin) / 2, 1) > 0
    error('There is something wrong with your parameter list. Probably missing a parameter value.');
end
count = zeros(length(pnames), 1); % variable to check if parameters are specified only once
    
for k=1:2:length(varargin)
    param = tolower(varargin{k});
    
    t = ismember(pnames,param);
    count = count + t;
    if sum(count > 1) > 0
        error(['Parameter "', param, '" specified more than once in the parameter list!']);
    end
    
    if sum(t) == 1
        params = setfield(params, param, varargin{k + 1});
    elseif sum(t) > 1
        error(['Parameter "', param, '" specified more than once in the default parameter list!']);
    else
        error(['Unknown parameter: "', param, '"']);
    end
end


%-------------------------
% Checking field types
%-------------------------
if isstruct(B)
    fn          = fieldnames(B);
    fname       = B.(fn{1});
    pth         = B.(fn{2});
    delfile     = 0;
elseif ischar(B)
    fname       = B;
    pth         = [];
    delfile     = 0;
elseif strcmp(params.ivartype,'cs')
    params.ivartype    = 'diag';
    [rows,cols] = size(B);
    if isequal(rows,cols) && isequal(rows,lmax+1)
        B = cssc2clm(B,lmax);
    elseif isequal(cols,(2*lmax+1)) && isequal(rows,lmax+1)
        B = cssc2clm(B,lmax)
    elseif ~isequal(cols,4) || ~isequal(rows, sum(1:lmax+1))
        error('Error matrix not in the SC, CS, or look-up table formats')
    end
    vec = sparse(diag([B(:,3); B(lmax+2:end,4)])); save 'qmatnow.mat' vec
    fname       = 'qmatnow.mat';
    pth         = pwd;
    pth         = [pth,'/'];
    delfile     = 1;
else
    [rows,cols] = size(B);
    if ~isequal(rows,cols) || ~isequal((lmax+1)^2,rows)
        error('Matrix dimensions do not agree with the maximum degree of the spherical harmonic expansion')
    end
    save 'qmatnow.mat' B
    fname       = 'qmatnow.mat';
    pth         = pwd;
    pth         = [pth,'/'];
    delfile     = 1;
end

%----------------------------------------------------
% Checking required input and output matrix types
%----------------------------------------------------
if strcmp(params.ivartype,'full')
    ivartype = 0;
elseif strcmp(params.ivartype,'block')
    ivartype = 1;
elseif strcmp(params.ivartype,'diag')
    ivartype = 2;
else
    display('Warning: Unknown covariance matrix structure')
end

fprintf('done\t %g[s] \n',toc)

tic
% ------------------------
% Grid definition.
% ------------------------
fprintf('Preparing for the desired grid ... ')
grd = params.grid;
loc = params.location;

if ~strcmp(grd,'points')
    if isempty(params.gridsize)
        dt = pi/180;
    elseif (params.gridsize < 0) || (params.gridsize > pi)
        error('Grid size value must be positive')
    else
        dt = params.gridsize;
    end
    
    if ~ischar(params.gridext)
        if length(params.gridext) == 1
            len = ones(4,1)*params.gridext;
        elseif length(params.gridext) == 2
            len = (params.gridext(:))';
            len = [params.gridext params.gridext];
        elseif length(params.gridext) == 3
            len = (params.gridext(:))';
            len = [params.gridext params.gridext(3)];
        end
        syn  = 'regional';
        % ext = [loc(:,1)-len(1), loc(:,1)+len(3), loc(:,2)-len(2), loc(:,2)+len(4), ones(size(loc))*dt];
        ext  = [loc(:,1)-len(1), loc(:,1)+len(3), loc(:,2)-len(2), loc(:,2)+len(4)];
        if strcmp(grd,'neumann') || strcmp(grd,'gauss')
            error('For regional computations Gauss(-Neumann) grid is not supported')
        end
    elseif strcmp(params.gridext,'global')
        syn = 'global';
    end
    
elseif strcmp(grd,'points')
    if size(params.gridext,2)~=2
        error('Parameter ''gridext'' does not contain points for the ''points'' grid type chosen')
    else
        syn = 'points';
        phi = pi/2 - params.gridext(:,1);
        lam = params.gridext(:,2);
        r   = 6378137;
    end
else
    error('Unknown grid type')
end
fprintf(' done \n')


%-----------------------------------
% Processing filter (coefficients)
%-----------------------------------
fprintf('Preparing the filter coefficients ...')

if ~ischar(params.filter)
    [rows,cols] = size(params.filter);
    
    if isequal(cols,2*lmax+1) && isequal(rows,lmax+1)
        filcff = css2clm(params.filter,lmax);
    elseif isequal(cols,4) && isequal(rows,sum(1:lmax+1))
        filcff = params.filter;
    elseif (cols==2) && (rows==lmax+1)
        filcff = cssc2clm(params.filter, lmax);
    elseif isequal(rows,cols) && isequal(rows,lmax+1)
        filcff = css2clm(params.filter,lmax);
    else
        error('The filter coefficients is not in the required format')
    end
else
    filcff = cssc2clm([(0:lmax)' ones(lmax+1,1)],lmax);
end

fprintf(' done \t %g[s]\n',toc)

tic
%---------------------------------------------------------------------
% Step 1: Compute spherical harmonics for the specified locations 'y'
%---------------------------------------------------------------------
fprintf('Computing spherical harmonics for the locations of the two-point functions ...')
[ylmc, ylms] = ylmdsgn(loc,lmax,params.quant,params.const);

% Treating the base functions with the filter coefficients
ylmc = bsxfun(@times, ylmc, filcff(:,3)');
ylms = [zeros(size(ylms,1),lmax+1) ylms];
ylms = bsxfun(@times, ylms, filcff(:,4)');
fprintf(' done\n')

fprintf('Starting spherical harmonic synthesis ')
% Generating indices for the beginning of the orders
idx2 = cumsum(lmax+1:-1:1);
idx1 = [1,idx2(1:end-1)+1];

Bc = zeros(size(ylmc));
Bs = Bc;

b = cell(size(loc,1),1);

cnt = unique(round(linspace(1,lmax+1,10)));

%-----------------------------------
% Step 2: Compute the product 'y*B'
%-----------------------------------
if ivartype==0

	fprintf('of fully populated matrix\n')
    fprintf('Reading fully populated matrix ... ')
    covcell = covord(fname,pth,lmax,0);
    if delfile == 1, delete(fname), end
    fprintf('done \t %g[s]\n',toc)
    
    tic
    fprintf('Computing y*B of y*B*Y'' ');

    for m = 1:lmax+1
        for k = m:lmax+1
            Bc(:,idx1(k):idx2(k)) = Bc(:,idx1(k):idx2(k)) + ...
                                    ylmc(:,idx1(m):idx2(m))*covcell.cc{m,k} + ...
                                    ylms(:,idx1(m):idx2(m))*covcell.sc{m,k};
            Bs(:,idx1(k):idx2(k)) = Bs(:,idx1(k):idx2(k)) + ...
                                    ylmc(:,idx1(m):idx2(m))*covcell.cs{m,k} + ...
                                    ylms(:,idx1(m):idx2(m))*covcell.ss{m,k};

            if m~=k
                Bc(:,idx1(m):idx2(m)) = Bc(:,idx1(m):idx2(m)) + ...
                                        ylmc(:,idx1(k):idx2(k))*covcell.cc{k,m} + ...
                                        ylms(:,idx1(k):idx2(k))*covcell.sc{k,m};
                Bs(:,idx1(m):idx2(m)) = Bs(:,idx1(m):idx2(m)) + ...
                                        ylmc(:,idx1(k):idx2(k))*covcell.cs{k,m} + ...
                                        ylms(:,idx1(k):idx2(k))*covcell.ss{k,m};
            end
        end
        if ismember(m,cnt)
            fprintf('.')
        end
    end
    fprintf(' done \t %g[s]\n',toc)

elseif ivartype==1
    
    fprintf('of block-diagonal matrix\n')
    fprintf('Reading block-diagonal matrix ... ')
    covcell = covord(fname,pth,lmax,1);
    if delfile == 1, delete(fname), end
    fprintf('done\t %g[s]\n',toc)
 
    tic
    fprintf('Computing y*B of y*B*Y'' ');

    for m = 1:lmax+1
        Bc(:,idx1(m):idx2(m)) = ylmc(:,idx1(m):idx2(m))*covcell.cc{m};
        Bs(:,idx1(m):idx2(m)) = ylms(:,idx1(m):idx2(m))*covcell.ss{m};

        if ismember(m,cnt)
            fprintf('.')
        end
    end
    fprintf(' done \t %g[s]\n',toc)


elseif ivartype == 2
    
    fprintf('of diagonal matrix\n')
    fprintf('Reading diagonal matrix ... ')
    covcell = covord(fname,pth,lmax,2);
    if delfile == 1, delete(fname), end
    fprintf('done\t %g[s] \n',toc)

    tic
    fprintf('Computing y*B of y*B*Y ');

    for m = 1:lmax+1
        Bc(:,idx1(m):idx2(m)) = bsxfun(@times,ylmc(:,idx1(m):idx2(m)),covcell.cc{m}');
        Bs(:,idx1(m):idx2(m)) = bsxfun(@times,ylms(:,idx1(m):idx2(m)),covcell.ss{m}');

        if ismember(m,cnt)
            fprintf('.')
        end
    end
    fprintf(' done \t %g[s]\n',toc)


end

tic
%----------------------------------------------------------------------------
% Step 3: Perform spherical harmonic synthesis over each of the points yB*Y'
%----------------------------------------------------------------------------
fprintf('Computing yB * Y'' of y*B*Y'' for \n')
Bc = bsxfun(@times,Bc,filcff(:,3)');
Bs = bsxfun(@times,Bs,filcff(:,4)');

switch params.gridcoord
case{'geo'}
    for p = 1:length(b)
        strng = ['Co-latitude ', num2str(loc(p,1)*180/pi),', Longitude ', num2str(loc(p,2)*180/pi)];
        fprintf('%s \n', strng)
        Bcs  = clm2sc([filcff(:,1:2) Bc(p,:)' Bs(p,:)']);

        switch syn
        case{'points'}
            b{p} = gshs_ptw(Bcs, lam, phi, r, 6378137, 'quant', params.quant, ...
                                                       'sub_WGS84', false, ...
                                                       'legendre', 'mex');
        case{'global'}
           [b{p}, th, lam] = gshs_(Bcs, 'quant', params.quant, ...
                                         'grid', grd, ...
                                         'gridsize', round(pi/dt), ...
                                         'height', params.height, ...
                                         'sub_WGS84', false, ...
                                         'legendre', 'mex');
        case{'regional'}
            % [b{p},th,lam] = shsreggrid(Bcs,ext(p,:),params.quant,grd,params.height,0);
            [b{p}, th, lam] = regionalshs(Bcs, ext(p,:), 'quant', params.quant, ...
                                                         'grid', grd, ...
                                                         'dlat', dt, ...
                                                         'dlong', dt, ...
                                                         'height', params.height, ...
                                                         'sub_normal', false, ...
                                                         'legendre', 'mex');
        end
    end

case{'polar'}
    for p = 1:length(b)
        fprintf('Co-latitude %g, Longitude %g \n', loc(p,1)*180/pi, loc(p,2)*180/pi)
        Bcs = clm2sc([filcff(:,1:2) Bc(p,:)' Bs(p,:)']);
        Bcs = real2cpxsh(Bcs,lmax);

        for k = 1:lmax
            C   = cpx2realmat(k);
            D   = bsxfun(@times,exp(1i*(-k:k)'*loc(p,2)),dlmk(k,-loc(p,1)));
            idx = lmax+1-k:lmax+1+k;
            Bcs(k+1,idx) = Bcs(k+1,idx)*D;
        end
        Bcs = cpx2realsh(Bcs,lmax);

        switch syn
        case{'points'}
            b{p} = gshs_ptw(Bcs, lam, phi, r, 6378137, 'quant', params.quant, ...
                                                       'sub_WGS84', false, ...
                                                       'legendre', 'mex');
        case{'global'}
            [b{p}, th, lam] = gshs_(Bcs, 'quant', params.quant, ...
                                         'grid', grd, ...
                                         'gridsize', round(pi/dt), ...
                                         'height', params.height, ...
                                         'sub_WGS84', false, ...
                                         'legendre', 'mex');
            lam = pi - lam;
        case{'regional'}
            row = size(loc,1);
            txe = [0 max(ext(p,2:3)) 0 2*pi];
            % [b{p},th,lam] = shsreggrid(Bcs,txe,params.quant,grd,params.height,0);
            [b{p}, th, lam] = regionalshs(Bcs, txe, 'quant', params.quant, ...
                                                    'grid', grd, ...
                                                    'dlat', dt, ...
                                                    'dlong', dt, ...
                                                    'height', params.height, ...
                                                    'sub_normal', false, ...
                                                    'legendre', 'mex');
            lam = pi - lam;
        end
    end
end
toc

fprintf('Synthesis of the two-point function complete. \n\n')
