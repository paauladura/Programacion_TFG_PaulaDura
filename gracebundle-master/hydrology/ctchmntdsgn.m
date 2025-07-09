function [A,indx] = ctchmntdsgn(cindx,lmax,grid,awt,quant,fltr,h,strct)

% CTCHMNTDSGN constructs a design matrix of aggregated catchment pixels for
% a given spherical harmonic expansion.
%
% A         = ctchmntdsgn(cindx,lmax)
% [A,indx]  = ctchmntdsgn(cindx,lmax,grid,awt,quant,fltr,h,strct)
% 
% INPUT
% cindx -   Catchment index grid. The grid must be from -180 to 180 in 
%           the longitude direction and 90 to -90 in latitude direction.
% lmax  -   Maximum expansion of the spherical harmonic coefficients
% grid  -   String argument, defining the catchment index grid:
%               1. 'pole' or 'mesh': equi-angular (N+1)*2N, including poles
%               and Greenwich meridian.                     -def: 'pole'
%               2. 'block' or 'cell': equi-angular block mid-points.
%               3. 'neumann' or 'gauss': Gauss-grid (N+1)*2N
% awt   -   Area-weight switch. If '1' then an area-weighted A matrix is
%           output. 
%                                                           -def: 0
% quant -   optional string argument, defining the field quantity:
%           'none' (coefficients define the output), 'geoid',
%           'potential', 'dg' or 'gravity' (grav. anomaly), 'tr' (grav.
%           disturbance), or 'trr' (2nd rad. derivative), 'water'
%           (water equivalent heights), 'smd' (surface mass density).
%                                                           -def: 'none'
% fltr  -   A structure containing the filter name and the filter parameters
%           It can handle only homogeneous isotropic filters.
%           'name' - A string with the one of the following filter names
%               1. Gauss
%               2. spCosine
%               3. Boxcar
%               4. Butterworth
%               5. Pellinen
%               6. Diffusion
%               7. Cosine (spectral analog of vonHann filter)
%               8. spButter (spatial analog of spectral Butterworth filter)
%               9. other (any other filter whose spectrum is known)
%                   Filter co-efficients should either be 
%                   given as a [l Wl] vector (isotropic case) or as a 
%                   matrix in CS- / SC- / [l m Clm Slm] in colombo 
%                   ordering formats (anisotropic case).
%
%           filter parameters for the filter provided in 'name'.
%           e.g. Gauss 500 km
%           name = 'Gauss'
%           cap = 500/ae, ae - semi-major axis of the reference ellipsoid
%
%           e.g. Spatial Cosine 200 km
%           name = 'spCosine'
%           cap = 200/ae
%
%           e.g. Box-car degree 60
%           name = 'Boxcar'
%           lc   = 60
%
%           e.g. Butterworth with cut-off degree 25 and order 5
%           name = 'Butterworth'
%           lc   = 25 (cut-off degree)
%           k    = 5 (order)
%           lmax = 200
%
%           e.g. Pellinen 800 km
%           name = 'Pellinen'
%           cap = 800/ae
%           lmax = 800
%
%           e.g. Diffusion with cut-off degree 40 and order 2
%           name = 'Diffusion'
%           lc   = 40 (cut-off degree)
%           k    = 2 (order)
%           lmax = 200
%
%           e.g. Spectral Cosine with start degree 8 cut-off degree 32 
%           and filter order 4
%           name    = 'Cosine'
%           ls      = 8
%           lc      = 32 (cut-off degree)
%           k       = 4 (order)
%           lmax    = 60
%
%           e.g. Spatial Butterworth with a cap of 1000 km and order 5
%           name = 'spButter'
%           cap  = 1000/ae
%           k    = 5 (order)
%
%           e.g. Any filter whose spectrum is known
%           name    = 'other'
%           spk     = Bl (Legendre spectrum of the filter)
%           lmax    = 90
%                                                           - def: 'none'
%
% h     -   height of evaluation of the field quantities in [m]
%                                                           -def: 0
% strct -   0. provides the design matrix
%           1. provides the design matrix in a structure variable
%           containing separate fields for cosine and sine components.
%           A = struct(A.c,A.s)
%           2. The cosine and sine fields have cells for each order of the
%           spherical harmonic coefficient expansion
%           A = struct('c',cell(1,lamx+1),'s',cell(1,lmax+1))
%                                                           -def: 0
%
% OUTPUT
% A     -   Design matrix provided in order-leading format
% indx  -   Unique index numbers of the catchments
%--------------------------------------------------------------------------
%
% See also blddsgn, ctchmntcov
%--------------------------------------------------------------------------

% USES CS2SC, CSSC2CLM, ISOTF

% Created on: 3 August 2009, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

if nargin < 3, error('Insufficient input arguments'),   end
if nargin == 8
    if isempty(strct),   strct = 0;          end
    if isempty(h),       h = 0;              end
    if isempty(fltr),    fltr = struct('type','none'); end
    if isempty(quant),   quant = 'none';     end
    if isempty(awt),     awt = 0;            end
elseif nargin == 7
    strct = 0;
    if isempty(h),       h = 0;              end
    if isempty(fltr),    fltr = struct('type','none'); end
    if isempty(quant),   quant = 'none';     end
    if isempty(awt),     awt = 0;            end
elseif nargin == 6
    strct = 0;          
    h = 0;
    if isempty(fltr),    fltr = struct('type','none'); end
    if isempty(quant),   quant = 'none';     end
    if isempty(awt),     awt = 0;            end
elseif nargin == 5
    strct = 0;          
    h = 0;              
    fltr = struct('type','none');
    if isempty(quant),   quant = 'none';     end
    if isempty(awt),     awt = 0;            end
elseif nargin == 4
    strct = 0;
    h = 0;
    fltr = struct('type','none');
    quant = 'none';
    if isempty(awt),     awt = 0;            end
elseif nargin == 3
    strct = 0; 
    h = 0;
    fltr = struct('type','none');
    quant = 'none';
    awt = 0;
end



%-------------------
% Useful constants
%-------------------
GM = 398600441800000;
ae = 6378136.6; % [m]

% ------------------------
% Grid definition.
% ------------------------
[n,m] = size(cindx);
cindx = [cindx(:,(floor(m/2) + 1):end) cindx(:,1:floor(m/2))];
if strcmp(grid,'pole') || strcmp(grid,'mesh')
    dt = 180/n-1;
    dl = 360/m;
    theta = (0:dt:180)';
    lam   = (0:dl:360-dt);
elseif strcmp(grid,'block') || strcmp(grid,'cell')
    dt = 180/n;
    dl = 360/m;
    theta = (dt/2:dt:180)';
    lam   = (dt/2:dl:360);
elseif strcmp(grid,'neumann') || strcmp(grid,'gauss')
    dt = 180/n-1;
    dl = 360/m;
    gx = grule(n);
    theta = flipud(acos(standing(gx)))*180/pi;
    lam   = (0:dl:360-dt);
else
    error('Unknown grid type')
end

%---------------------------------------------
% Calculating spectral filter coefficients
%---------------------------------------------

if isstruct(fltr)
    fldnm = fieldnames(fltr);
    
    switch fltr.name
        case {'Gauss', 'gauss', 'Gaussian', 'gaussian'}
            if ismember('cap',fldnm)
                B = gaussfltr(fltr.cap,lmax);
            else
                B = gaussfltr(pi/36,lmax);
            end
            B = B(1:lmax+1);
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case {'vonHann', 'vonhann', 'hann', 'Hann'}
            if ismember('cap',fldnm)
                B = vonhann(fltr.cap,lmax);
            else
                B = vonhann(pi/18,lmax);
            end
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case {'spCosine', 'spcosine', 'spCos', 'spcos'}
            if ismember('cap',fldnm)
                cap    = fltr.cap;
            else
                cap    = pi/18;
            end
            if ismember('k',fldnm)
                k      = fltr.k;
            else
                k      = 2;
            end
            B = sptcosine(cap,k,'Gauss',lmax,true);
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case 'Pellinen'
            perf.fil = 'Pellinen';
            if any(strcmp('cap',fldnm))
                cap    = fltr.cap;
            else
                cap    = pi/18;
            end
            B = pellinen(cap,(0:lmax)');
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case 'Box-car'
            if any(strcmp('lc',fldnm))
                lc = fltr.lc;
            else
                lc = 20;
            end
            B = [ones(lc+1,1); zeros(lmax-lc,1)];
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case 'Butterworth'
            if any(strcmp('lc',fldnm))
                lc = fltr.lc;
            else
                lc = 20;
            end
            if any(strcmp('k',fldnm))
                k      = fltr.k;
            else
                k      = 2;
            end
            B = bttrwrth(lc,k,lmax);
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case 'Diffusion'
            if any(strcmp('lc',fldnm))
                lc = fltr.lc;
            else
                lc = 20;
            end
            if any(strcmp('k',fldnm))
                k      = fltr.k;
            else
                k      = 1;
            end
            B = diffusionfil(lc,k,lmax);
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case {'Cosine', 'cosine', 'cos', 'Cos', 'spkcos', 'spkCos'}
            if any(strcmp('lc',fldnm))
                lc = fltr.lc;
            else
                if lmax > 60
                    lc = 60;
                else
                    lc = lmax;
                end
            end
            if any(strcmp('k',fldnm))
                k      = fltr.k;
            else
                k      = 1;
            end
            if any(strcmp('k',fldnm))
                ls      = fltr.ls;
            else
                ls      = 0;
            end
            B = spkcosine(lc,k,ls,lmax);
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case {'spButter', 'spbutter', 'Spbutter', 'SpButter'}
            if any(strcmp('cap',fldnm))
                cap    = fltr.cap;
            else
                cap    = pi/18;
            end
            if any(strcmp('k',fldnm))
                k      = fltr.k;
            else
                k      = 2;
            end
            B = sptbttrwrth(cap,k,'Gauss',lmax,true);
            fltrcff = cssc2clm([(0:lmax)' B],lmax);
        case 'other'
            if ~any(strcmp('spk',fldnm))
                error('Insufficient input arguments')
            else
                B = fltr.spk;
            end
            [r_,c_] = size(B);
            
            if isequal(c_,2*lmax+1) && isequal(r_,lmax+1)
                filcff = css2clm(B,lmax);
            elseif isequal(c_,4) && isequal(r_,sum(1:lmax+1))
                filcff = B;
            elseif (c_==2) && (r_==lmax+1)
                filcff = cssc2clm(B, lmax);
            elseif isequal(r_,c_) && isequal(r_,lmax+1)
                filcff = css2clm(B,lmax);
            else
                error('The filter coefficients is not in the required format')
            end
        otherwise
            error('String not recognized. Please check the filter name given.')
    end
    % fldvar = fieldnames(fltr);
    % if strcmp('gauss',fltr.(fldvar{1}))
    %     cap     = fltr.(fldvar{2});
    %     fltrcff = cssc2clm([(0:lmax)' ones(lmax+1,1)],lmax);
    % elseif strcmp('pell',fltr.(fldvar{1}))
    %     cap     = fltr.(fldvar{2});
    %     fltrcff = cssc2clm([(0:lmax)' ones(lmax+1,1)],lmax);
    % elseif strcmp('han',fltr.(fldvar{1}))
    %     cap     = 0;
    %     fltrcff = hannoniso(lmax,fltr.(fldvar{2}),fltr.(fldvar{3}),fltr.(fldvar{4}));
    % elseif strcmp('hann',fltr.(fldvar{1}))
    %     cap     = 0;
    %     fltrcff = vonhann(lmax,fltr.(fldvar{2}));
    %     fltrcff = cssc2clm([(0:lmax)' fltrcff],lmax);
    % elseif strcmp('other',fltr.(fldvar{1}))
    %     cap     = 0;
    %     fltrcff = cssc2clm(fltr.(fldvar{2}),lmax);
    % elseif strcmp('none',fltr.(fldvar{1}))
    %     fltrcff = cssc2clm([(0:lmax)' ones(lmax+1,1)],lmax);
    %     cap     = 0;
    % else
    %     error('Unknown filter type')
    % end
elseif ischar(fltr)
    if strcmp('none',fltr)
        fltrcff = cssc2clm(ones(lmax+1),lmax);
        cap     = 0;
    else
        error('Unknown filter type')
    end
end

%------------------------------------------------
% Calculating isotropic transfer coefficients
%------------------------------------------------
transf = eigengrav((0:lmax)', quant, h, [GM, ae]);
transf = cssc2clm([(0:lmax)' transf], lmax);
transf(:,3) = transf(:,3).*fltrcff(:,3);
transf(:,4) = transf(:,4).*fltrcff(:,4);

indx = unique(cindx(cindx(:) ~= -9999));

A = struct('c',zeros(length(indx),sum(1:lmax+1)),'s',zeros(length(indx),sum(1:lmax+1)));

if awt == 1
    ar = ae^2 * dl * dt * sind(theta');
elseif awt == 0
    ar = ones(size(theta'));
else
    error('Unknown switch value for area weights')
end

for i = 1:length(indx)
    [j,k] = find(cindx==indx(i));
    c = (transf(:,2)*lam(k))';
    s = sind(c);
    c = cosd(c);
    p = zeros(length(j),sum(1:lmax+1));
    mlen2 = cumsum(fliplr(1:lmax+1))';
    mlen1 = [1;mlen2(1:end-1)+1];
    for cnt = 1:lmax+1
        p(:,mlen1(cnt):mlen2(cnt)) = plm(cnt-1:lmax,cnt-1,theta(j)*pi/180);
    end
    w = ar(j)/sum(ar(j));
    A.c(i,:) = (w*(p.*c)).*transf(:,3)';
    A.s(i,:) = (w*(p.*s)).*transf(:,4)';
end

if strct == 0,
    A = [A.c A.s(:,(lmax+2):end)];
elseif strct == 1
    A.s = A.s(:,(lmax+2):end);
elseif strct == 2
    A.c = mat2cell(A.c,length(indx),(lmax+1:-1:1));
    A.s = mat2cell(A.s,length(indx),(lmax+1:-1:1));
end
