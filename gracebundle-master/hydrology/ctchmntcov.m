function [f, indx] = ctchmntcov(cindx, field, lmax, ivartype, quant, grd, fil, h)

% CTCHMNTCOV propagates the covariance information of a set of spherical 
% harmonic coefficients to the catchments specified by catchment index file
% by area-weighted aggregation. Error propagation is possible for all types
% of structures of the covariance matrix of the spherical harmonic
% coefficients.
%
% fvar        = ctchmntcov(cindx,field,lmax,ivartype)
% [fvar,indx] = ctchmntcov(cindx,field,lmax,ivartype,quant,grd,fil,h)
%
% INPUT
% cindx     -   Catchment index grid. The grid must be from -180 to 180 in 
%               the longitude direction and 90 to -90 along the latitude.
% field     -   A [(L+1)^2 * (L+1)^2] variance-covariance matrix arranged
%               in Colombo ordering, or, a structure with two variables:
%               1. a character array with the covariance matrix filename;
%               2. path where the file is stored;
%               or, a character array of the filename and the path
%               together,
%               or, if only the variances are available then provide them
%               in either SC, CS, or Colombo ordering formats
%               Use the first option only when the maximum degree of the
%               spherical harmonic expansion is [lmax <= 70].
%               Only the lower triangular part of the matrix must be used.
% lmax      -   Maximum degree of the spherical harmonic expansion.
% quant     -   optional string argument, defining the field quantity:
%               'none' (coefficients define the output), 'geoid',
%               'potential', 'dg' or 'gravity' (grav. anomaly), 'tr' (grav.
%               disturbance), or 'trr' (2nd rad. derivative), 'water'
%               (water equivalent heights), 'smd' (surface mass density).
%               If in the case of cross-covariance computation quant is a
%               cell array with two-variables. One for each field.
%                                                           -def: 'none'
% grd       -   optional string argument, defining the grid:
%               1. 'pole' or 'mesh': equi-angular (N+1)*2N, including poles
%               and Greenwich meridian.                     -def: 'pole'
%               2. 'block' or 'cell': equi-angular block mid-points.
%               3. 'neumann' or 'gauss': Gauss-grid (N+1)*2N
% fil       -   A structure containing the filter name and the filter parameters
%               It can handle only homogeneous isotropic filters.
%               'name' - A string with the one of the following filter names
%                   1. Gauss
%                   2. spCosine
%                   3. Boxcar
%                   4. Butterworth
%                   5. Pellinen
%                   6. Diffusion
%                   7. Cosine (spectral analog of vonHann filter)
%                   8. spButter (spatial analog of spectral Butterworth filter)
%                   9. other (any other filter whose spectrum is known)
%                       Filter co-efficients should either be 
%                       given as a [l Wl] vector (isotropic case) or as a 
%                       matrix in CS- / SC- / [l m Clm Slm] in colombo 
%                       ordering formats (anisotropic case).
%
%               filter parameters for the filter provided in 'name'.
%               e.g. Gauss 500 km
%               name = 'Gauss'
%               cap = 500/ae, ae - semi-major axis of the reference ellipsoid
%
%               e.g. Spatial Cosine 200 km
%               name = 'spCosine'
%               cap = 200/ae
%
%               e.g. Box-car degree 60
%               name = 'Boxcar'
%               lc   = 60
%
%               e.g. Butterworth with cut-off degree 25 and order 5
%               name = 'Butterworth'
%               lc   = 25 (cut-off degree)
%               k    = 5 (order)
%               lmax = 200
%
%               e.g. Pellinen 800 km
%               name = 'Pellinen'
%               cap = 800/ae
%               lmax = 800
%
%               e.g. Diffusion with cut-off degree 40 and order 2
%               name = 'Diffusion'
%               lc   = 40 (cut-off degree)
%               k    = 2 (order)
%               lmax = 200
%
%               e.g. Spectral Cosine with start degree 8 cut-off degree 32 
%               and filter order 4
%               name    = 'Cosine'
%               ls      = 8
%               lc      = 32 (cut-off degree)
%               k       = 4 (order)
%               lmax    = 60
%
%               e.g. Spatial Butterworth with a cap of 1000 km and order 5
%               name = 'spButter'
%               cap  = 1000/ae
%               k    = 5 (order)
%
%               e.g. Any filter whose spectrum is known
%               name    = 'other'
%               spk     = Bl (Legendre spectrum of the filter)
%               lmax    = 90
%                                                           - def: 'none'
% h         -   height for upward continuation [m]. It is also possible to
%               provide to different heights if the covariance matrix is
%               between two fields that are representative of two different
%               heights.
%                                                           - def: 0.
% ivartype  -   type of variance-covariance propagation
%               1. 'full' - full spectral variance-covariance matrix is
%                           propagated.
%               2. 'block'- block-diagonal variance-covariance matrix is
%                           propagated.
%               3. 'diag' - only the diagonal of the variance-covariance
%                           matrix is propagated.
%               4. 'cssc' - if only the diagonal elements are available and
%                           provided in CS/SC formats
%
% OUTPUT
% fvar      -   the global variances of the field
% indx      -   Catchment index vector
%--------------------------------------------------------------------------
%
% See also gshscov gshscovfn ctchmntdsgn
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Uses CS2SC, CTCHMNTDSGN, COVORD
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Created on: 31 July 2009, Stuttgart
% Author: Balaji Devaraju
%------------------------------------------------------

tic
display('-----------------------------------------------------------------')
display('      Error covariance propagation to catchment time-series')
display('-----------------------------------------------------------------')

%-----------------------------------------------------
% Checking the input and augmenting initial values
%-----------------------------------------------------
if nargin > 8, error('Too many input arguments!'),          end
if nargin == 7 || isempty(h),        h        = 0;             end
if nargin == 6 || isempty(fil),      fil      = struct('type','none'); end
if nargin == 5 || isempty(grd),      grd     = 'mesh';       end
if nargin == 4 || isempty(quant),    quant    = 'none';       end
if nargin == 3, error('Insufficient input arguments'),       end

if ~ischar(grd), error('grid argument must be string'),    end
grd  = lower(grd);

%-------------------------
% Checking field types
%-------------------------
if isstruct(field)
    fn          = fieldnames(field);
    fname       = field.(fn{1});
    pth         = field.(fn{2});
    delfile     = 0;
elseif ischar(field)
    fname       = field;
    pth         = [];
    delfile     = 0;
elseif strcmp(ivartype,'cssc')
    ivartype = 'diag';
    vec = cssc2clm(field,lmax);
    vec = sparse(diag([vec(:,3); vec(lmax+2:end,4)]));
    save 'qmatnow.mat' vec
    vec         = [];
    fname       = 'qmatnow.mat';
    pth         = pwd;
    pth         = [pth,'/'];
    delfile     = 1;
else
    [rows,cols] = size(field);
    if ~isequal(rows,cols) || ~isequal((lmax+1)^2,rows)
        error('Matrix dimensions do not agree with the maximum degree')
    end
    save 'qmatnow.mat' field
    fname       = 'qmatnow.mat';
    pth         = pwd;
    pth         = [pth,'/'];
    delfile     = 1;
end

%----------------------------------------------------
% Checking required input and output matrix types
%----------------------------------------------------
if strcmp(ivartype,'full')
    ivartype = 0;
elseif strcmp(ivartype,'block')
    ivartype = 1;
elseif strcmp(ivartype,'diag')
    ivartype = 2;
else
    display('Warning: Unknown covariance matrix structure')
end

toc

display('Checking of input data complete. Starting preliminary steps')

tic

%-----------------------------------------
% Calculating the catchment design matrix
%-----------------------------------------
[A,indx] = ctchmntdsgn(cindx,lmax,grd,1,quant,fil,h,2);


toc

display('Preliminary steps completed. Starting Spherical harmonic synthesis.')

tic

hwb = twaitbar('init', [], 'gshs');  % initialize the waitbar
if ivartype==0    
    covcell = covord(fname,pth,lmax,0);
    if delfile == 1, delete(fname), end
    f = 0;
    
    for k = 1:lmax+1
        for m = k:lmax+1
            if m==k
                covcell.cc{m,k} = covcell.cc{m,k} + covcell.cc{m,k}' - diag(diag(covcell.cc{m,k}));
                covcell.ss{m,k} = covcell.ss{m,k} + covcell.ss{m,k}' - diag(diag(covcell.ss{m,k}));
            end
            f = f + A.c{m}*covcell.cc{m,k}*A.c{k}';
            f = f + A.c{m}*covcell.sc{k,m}'*A.s{k}';
            f = f + A.s{m}*covcell.sc{m,k}*A.c{k}';
            f = f + A.s{m}*covcell.ss{m,k}*A.s{k}';
            if m~=k
                f = f + A.c{k}*covcell.cc{m,k}'*A.c{m}';
                f = f + A.c{k}*covcell.sc{m,k}'*A.s{m}';
                f = f + A.s{k}*covcell.sc{k,m}*A.c{m}';
                f = f + A.s{k}*covcell.ss{m,k}'*A.s{m}';
            end
        end
        hwb = twaitbar(k/(lmax + 1), hwb); % update the waitbar
    end
elseif ivartype==1
    covcell = covord(fname,pth,lmax,1);
    if delfile == 1, delete(fname), end
    f = 0;
    for m = 1:lmax+1
        covcell.cc{m} = covcell.cc{m} + covcell.cc{m}' - diag(diag(covcell.cc{m}));
        covcell.ss{m} = covcell.ss{m} + covcell.ss{m}' - diag(diag(covcell.ss{m}));
        f = f + A.c{m}*covcell.cc{m}*A.c{m}' + A.s{m}*covcell.ss{m}*A.s{m}';          
        hwb = twaitbar(m/(lmax+1), hwb);
    end
elseif ivartype == 2
    covcell = covord(fname,pth,lmax,2);
    if delfile == 1, delete(fname), end
    f = 0;
    for m = 1:lmax+1
        f = f + A.c{m}*diag(covcell.cc{m})*A.c{m}' + A.s{m}*diag(covcell.ss{m})*A.s{m}';          
        hwb = twaitbar(m/(lmax+1), hwb);
    end
end
twaitbar('close', hwb); % finalize the waitbar

toc
