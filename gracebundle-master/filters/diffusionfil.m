function W = diffusionfil(lc,k,lmax)

% DIFFUSIONFIL computes the spectral coefficients of the diffusion filter
% proposed by Sardeshmukh and Hoskins (1984), Spatial smoothing on the
% sphere, Monthly Weather Review, 112: 2524-2529.
%
% W = diffusionfil(lc,k,lmax)
%
% INPUT
% lc    -   Cut-off degree [positive integer]
% k     -   Tuning parameter [positive integer]
% lmax  -   Maximum degree of spherical harmonic expansion
%
% OUTPUT
% W     -   Spectral coefficients
%--------------------------------------------------------------------------

% Created on: 16 March 2010, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

if nargin == 0
    lc      = 40;
    k       = 2;
    lmax    = 120;
elseif nargin == 1
    if ~isint(lc) || (lc < 0)
        fprintf('WARNING: LC must be a real positive integer.')
        lc = fix(abs(lc));
    end
    k       = 2;
    lmax    = 2*lc;
elseif nargin == 2
    if ~isint(lc) || (lc < 0)
        fprintf('WARNING: LC must be a real positive integer.')
        lc  = fix(abs(lc));
    end
    if ~isint(k) || (k < 0)
        fprintf('WARNING: K must be a real positive integer.')
        k   = fix(abs(k));
    end
    lmax    = 2*lc;
elseif nargin == 3
    if ~isint(lc) || (lc < 0)
        fprintf('WARNING: LC must be a real positive integer.')
        lc  = fix(abs(lc));
    end
    if ~isint(k) || (k < 0)
        fprintf('WARNING: K must be a real positive integer.')
        k   = fix(abs(k));
    end
    if ~isint(lmax) || (k < 0)
        fprintf('WARNING: K must be a real positive integer.')
        lmax = fix(abs(lmax));
    end
end


W = 0:lmax;
W = exp(-((W.*(W+1))/(lc*(lc+1))).^k);
W = W(:);
