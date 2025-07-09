function varargout = sptcosine(psi0,k,psi,lmax,spk)

% SPTCOSINE computes the spatial coefficients of the cosine window of
% order k. The cosine window of order 2 is the von Hann window.
%
% w = cos(psi/psi0 * pi/2)^k, psi<=psi0; b = 0, psi>psi0
%
% [w,psi]   = sptcosine
% w/W       = sptcosine(psi0,k,psi,lmax)
% [W,w,psi] = sptcosine(psi0,k,psi,lmax,spk)
%
%
% INPUT
% psi0  - Filter radius in radians
% k     - Order of the spatial cosine filter [integer]
% psi   - Spherical distances at which the filter value is sought. If the
%         filter values are sought on Gauss-Neumann sampling points then
%         'Gauss'. The default is 181 sampling points. For a desired number
%         of sampling points provide lmax as well.
% lmax  - Maximum number of Legendre polynomial degree to generate the 
%         spectrum of the spatial cosine function. Following are the legal
%         values of LMAX: 
%           a) if LMAX is logical and true, generates spectrum with 
%               LMAX = fix(length(psi)/2), but if you specified 'Gauss' for 
%               PSI then LMAX = 180
%           b) if LMAX is an integer then the spectrum has LMAX degrees, only
%               in the case PSI is a vector of sampling points, else LMAX will
%               be used as the number of sampling points for the Gauss-Neumann
%               sampling.
% spk   - Spectrum estimation toggle switch. Use this only when you need 
%         Gauss-Neumann sampling and need to specify the number of samples
%         via LMAX. 
%
% OUTPUT
% w     - Spatial function values
% W     - Spectrum of spatial cosine function
% psi   - Sampling points [radians]
%
%
%
%--------------------------------------------------------------------------
% USES uberall/grule
%
% See also spkcosine, vonhann
%--------------------------------------------------------------------------
%

% Created on: 30 November 2007, Stuttgart
% Authors: Balaji Devaraju

% Revision History:
% BD    02/04/2014 Completely rewritten
%-------------------------------------------------------------------------------



if nargin == 0
    psi0    = pi/36;
    k       = 5;
    psi     = (0:pi/180:pi);
    spk     = false;
elseif nargin == 1
    k   = 5;
    psi = (0:pi/180:pi);
    spk = false;
elseif nargin == 2    
    psi = (0:pi/180:pi);
    spk = false;

    if isscalar(k)
        k = fix(k);
    else
        k = 2;
        fprintf('WARNING: variable "k" is not a scalar integer\n')
    end

    if ~isscalar(psi0)
        error('variable "psi0" must be a scalar')
    end
elseif nargin == 3
    if ischar(psi) && strcmp(psi,'Gauss')
        lmax        = 180;
        [tmp,wf]    = grule(lmax+1);
        wf          = flipud(wf(:));
        psi         = acos(flipud(tmp(:)));
    else
        psi = psi(:);
    end
    spk = false;

    if isscalar(k)
        k = fix(k);
    else
        k = 2;
        fprintf('WARNING: variable "k" is not a scalar integer\n')
    end

    if ~isscalar(psi0)
        error('variable "psi0" must be a scalar')
    end
elseif nargin == 4
    if islogical(lmax) && ~lmax
        spk     = false;
    elseif islogical(lmax) && lmax
        spk     = true;
    elseif ~islogical(lmax) && isscalar(lmax)
        lmax    = fix(lmax);
        spk     = false;
    end

    if ischar(psi) && strcmp(psi,'Gauss') 
        tmp = grule(lmax+1);
        psi = acos(flipud(tmp(:)));
        if spk
            lmax = 180;
        end
    elseif isvector(psi)
        psi = psi(:);
        if spk
            lmax = fix(length(psi/2));
        end
    elseif isscalar(psi)
        psi = fix(psi);
        psi = (0:pi/psi:pi)';
    end

    if isscalar(k)
        k = fix(k);
    else
        k = 2;
        fprintf('WARNING: variable "k" is not a scalar integer\n')
    end

    if ~isscalar(psi0)
        error('variable "psi0" must be a scalar')
    end
elseif nargin == 5 && islogical(spk)
    if ~islogical(spk)
        error('variable "spk" must be logical')
    end

    if ~isscalar(lmax)
        error('variable "lmax" must be a scalar')
    else
        lmax = fix(lmax);
    end

    if ischar(psi) && strcmp(psi,'Gauss')
        tmp = grule(lmax+1);
        psi = acos(flipud(tmp(:)));
    elseif isvector(psi)
        psi = psi(:);
    elseif isscalar(psi)
        psi = fix(psi);
        psi = (0:pi/psi:pi)';
    end

    if isscalar(k)
        k = fix(k);
    else
        k = 2;
        fprintf('WARNING: variable "k" is not a scalar integer\n')
    end

    if ~isscalar(psi0)
        error('variable "psi0" must be a scalar')
    end
else 
    error('Please verify your inputs.')
end

w 			    = (cos(psi/psi0 * pi/2)).^k;
w(psi > psi0) 	= 0;

if spk
    varargout{2}    = w;
    varargout{3}    = psi;

    if psi0 < (4*pi/lmax)
        L    = lmax;
        lmax = round(4*pi/psi0);
    end
 
    [tmp,wf]        = grule(lmax+1);
    tmp             = acos(flipud(tmp(:)));
    wf              = wf(:)/2;

    P               = legpol((0:lmax)',tmp);
    w               = (cos(tmp/psi0 * pi/2)).^k;
    w(tmp > psi0)   = 0;
    w               = w.*wf;
    W               = P'*w;

    if exist('L','var')
        W = W(1:L+1);
    end

    varargout{1}    = W/W(1);
else
    varargout{1}    = w;
    varargout{2}    = psi;
end
