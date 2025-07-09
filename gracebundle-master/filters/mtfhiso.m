function [mtf,b1,b2,th,sig1,sig2] = mtfhiso(fil,n,sep,start,fin)

% MTFHISO computes the Modulation Transfer Function for a homogeneous isotropic filter
%
% mtf = mtfhiso(fil,n,sep,start,fin)
%
% INPUT
% fil   -   A structure containing the filter name and the filter parameters
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
%
%
% n     - sampling points in [0,pi] co-latitudes
% sep   - The step size for increasing the separation between the dirac pulses
% start - Starting separation
% fin   - Separation length for which the MTF must be calculated
%
% OUTPUT
% mtf   - modulation transfer function [Sph.Distance ModTrans]
% b1,b2 - Filtered Dirac pulses: b1 --> [n x 1]
%                                b2 --> [n x (fin-start)/sep]
% th    - Sampling points
% sig1  - Location of the Dirac pulse 1
% sig2  - Locations of the Dirac pulse 2
%
%-------------------------------------------------------------------------
% USES FilterBundle uberall
%-------------------------------------------------------------------------

% Balaji Devaraju. Stuttgart, 10 June 2013.

if nargin < 1
    error('Insufficient input argumnts')
end

if ~isstruct(fil)
    error('The input fil must be a structure variable')
end

if nargin == 1
    n       = 7200;
    sep     = pi/1800;
    start   = pi/360;
    fin     = pi/4;
elseif nargin == 2
    sep     = pi/1800;
    start   = pi/360;
    fin     = pi/4;
elseif nargin == 3
    start   = pi/360;
    fin     = pi/4;
elseif nargin == 4
    fin     = pi/4;
end

[fbar,th,sig1,sig2,b1,b2] = fsmooth(fil,n,sep,start,fin);

mtf = zeros(length(sig2),2);
mtf(:,1) = sig2 - sig1;
for k = 1:size(fbar,2)
    mx = max(fbar(:,k));
    dt = abs(th-sig1-mtf(k,1)/2);
    ind = find(dt==min(dt));
    if ~isempty(ind)
        mtf(k,2) = 1 - fbar(ind(1),k)/mx;
    else
        error('There is a bug in finding the midpoint between the two Dirac pulses')
    end
end

end

%------------------------------------------------------------------------------
%                               FSMOOTH
%------------------------------------------------------------------------------

function [fbar,th,sig1,sig2,b1,b2] = fsmooth(fil,n,sep,start,fin)

% FSMOOTH provides the smoothed field to compute the MTF

% constants % loading all the constants

sig1 	= pi/6; % Dirac pulse location
sig2 	= sig1 + (start:sep:fin); % Dirac pulse location
th 		= (0:pi/n:(sig2(end)+pi/6))'; % smaple points

r = length(th);
c = length(sig2);

fldnm = fieldnames(fil);

sphd1   = abs(th - sig1);
sphd2   = abs(th*ones(1,c) - ones(r,1)*sig2);

switch lower(fil.name)
    case { 'gauss', 'gaussian' }
        if ~any(strcmp('cap',fldnm))
            fil.cap    = pi/36;
        end
        fil.cap = fil.cap;
        [b1,~]	= sptgauss(fil.cap,sphd1,1);
        [b2,~]	= sptgauss(fil.cap,sphd2(:),1);
        b2      = reshape(b2,r,c);
        
        fbar    = bsxfun(@plus,b1,b2);
        
    case { 'sptcosine', 'sptcos', 'spcos', 'spcosine', 'spatialcos' }
        if ~any(strcmp('cap',fldnm))
            fil.cap    = pi/18;
        end
        if isempty(fil.k)
            fil.k    = 2;
        end
        b1	= sptcosine(fil.cap,fil.k,sphd1);
        b2	= sptcosine(fil.cap,fil.k,sphd2(:));
        b2  = reshape(b2,r,c);
        
        fbar    = bsxfun(@plus,b1,b2);
        
    case { 'pellinen', 'pell' }
        if ~any(strcmp('cap',fldnm))
            fil.cap = pi/18;
        end
        if ~any(strcmp('lmax',fldnm))
            fil.lmax = 800;
        end
        % B       = pellinen(fil.cap,(0:fil.lmax)');
        % b1      = spk2spt(B,fil.lmax,sphd1);
        % fbar    = zeros(r,c);
        % for k = 1:c
        %     b2(:,k)     = spk2spt(B,fil.lmax,sphd2(:,k));
        %     fbar(:,k)   = b1 + b2(:,k);
        % end
        
        b1      = (sphd1<=fil.cap);
        b2      = (sphd2<=fil.cap);
        
        fbar    = bsxfun(@plus,b1,b2);
        
    case { 'box_car', 'box-car', 'boxcar', 'shannon' }
        if ~any(strcmp('lc',fldnm))
            fil.lc = 20;
        end
        b1	    = spk2spt(ones(fil.lc+1,1),fil.lc,sphd1);
        fbar    = zeros(r,c);
        b2      = fbar;
        for k = 1:c
            b2(:,k)     = spk2spt(ones(fil.lc+1,1),fil.lc,sphd2(:,k));
            fbar(:,k)   = b1 + b2(:,k);
        end
        
    case { 'butter', 'butterworth' }
        if ~any(strcmp('lc',fldnm))
            fil.lc = 20;
        end
        
        if ~any(strcmp('k',fldnm))
            fil.k = 5;
        end
        B       = bttrwrth(fil.lc,fil.k,fil.lmax);
        b1      = spk2spt(B,fil.lmax,sphd1);
        fbar    = zeros(r,c);
        b2      = fbar;
        for k = 1:c
            b2(:,k)     = spk2spt(B,fil.lmax,sphd2(:,k));
            fbar(:,k)   = b1 + b2(:,k);
        end
        
    case { 'diffusion' }
        if ~any(strcmp('lc',fldnm))
            fil.lc = 20;
        end
        if ~any(strcmp('k',fldnm))
            fil.k = 5;
        end
        if ~any(strcmp('lmax',fldnm))
            fil.lmax = 800;
        end
        B       = diffusionfil(fil.lc,fil.k,fil.lmax);
        b1      = spk2spt(B,fil.lmax,sphd1);
        fbar    = zeros(r,c);
        b2      = fbar;
        for k = 1:c
            b2(:,k)     = spk2spt(B,fil.lmax,sphd2(:,k));
            fbar(:,k)   = b1 + b2(:,k);
        end
        
    case { 'cosine', 'cos'}
        if ~any(strcmp('lc',fldnm))
            fil.lc = 20;
        end
        if ~any(strcmp('k',fldnm))
            fil.k = 5;
        end
        if ~any(strcmp('lmax',fldnm))
            fil.lmax = 800;
        end
        if ~any(strcmp('ls',fldnm))
            fil.ls = 0;
        end
        B       = spkcosine(fil.lc,fil.k,fil.ls,fil.lmax);
        b1      = spk2spt(B,fil.lmax,sphd1);
        fbar    = zeros(r,c);
        b2      = fbar;
        for k = 1:c
            b2(:,k)     = spk2spt(B,fil.lmax,sphd2(:,k));
            fbar(:,k)   = b1 + b2(:,k);
        end
        
    case { 'spbutter' }
        if ~any(strcmp('cap',fldnm))
            fil.cap = pi/18;
        end
        if ~any(strcmp('k',fldnm))
            fil.k = 5;
        end
        b1      = sptbttrwrth(fil.cap,sphd1,fil.k);
        b2      = sptbttrwrth(fil.cap,sphd2(:),fil.k);
        b2      = reshape(b2,r,c);
        fbar    = b1*ones(1,c) + b2;
        
    case 'other'
        if ~any(strcmp('spk',fldnm)) || ~any(strcmp('lmax',fldnm))
            error('Insufficient input arguments')
        else
            B       = fil.spk;
            b1      = spk2spt(B,fil.lmax,sphd1);
            fbar    = zeros(r,c);
            b2      = fbar;
            for k = 1:c
                b2(:,k)     = spk2spt(B,fil.lmax,sphd2(:,k));
                fbar(:,k)   = b1 + b2(:,k);
            end
        end
    otherwise
        error('Unknown filter type')
        
end
end
