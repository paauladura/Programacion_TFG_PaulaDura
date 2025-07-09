function varargout = sptbttrwrth(psi0,k,psi,lmax,spk)

% SPTBTTRWRTH generates the weights of the spatial Butterworth filter,
% which is the spatial analog of the spectral Butterworth filter on the
% sphere.
%
% [w,psi]   = sptbttrwrth
% w/W       = sptbttrwrth(psi0,k,psi,lmax)
% [W,w,psi] = sptbttrwrth(psi0,k,psi,lmax,spk)
%
% INPUT
% psi0  - Filter radius in radians
% k     - Order of the Butterworth filter [integer]
% psi   - Spherical distances at which the filter value is sought. If the
%         filter values are sought on Gauss-Neumann sampling points then
%         'Gauss'. The default is 181 sampling points. For a desired number
%         of sampling points provide lmax as well.
% lmax  - Maximum number of Legendre polynomial degree to generate the 
%         spectrum of the spatial Butterworth function. Following are legal
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
% W     - Spectrum of spatial Butterworth function
% psi   - Sampling points [radians]
%
%
%
%--------------------------------------------------------------------------
% USES uberall/grule
%
% See also bttrwrth
%--------------------------------------------------------------------------
%
%

% 10 May 2011, Stuttgart. Balaji Devaraju.

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
elseif nargin == 3
    if ischar(psi) && strcmp(psi,'Gauss')
        lmax        = 180;
        [tmp,wf]    = grule(lmax+1);
        wf          = flipud(wf(:));
        psi         = acos(flipud(tmp(:)));
        spk         = false;
    else
        psi = psi(:);
        spk = false;
    end
elseif nargin == 4
    if ischar(psi) && strcmp(psi,'Gauss')
        if islogical(lmax) && ~lmax
            spk     = false;
            lmax    = 180;
        elseif islogical(lmax) && lmax
            spk     = true;
            lmax    = 180;
        elseif ~islogical(lmax) && isscalar(lmax)
            lmax    = fix(lmax);
            spk     = false;
        end
        [tmp,wf]    = grule(lmax+1);
        wf          = flipud(wf(:))/2;
        psi         = acos(flipud(tmp(:)));
    elseif ~isempty(psi)
        if isvector(psi)
            psi = psi(:);
        elseif isscalar(psi)
            psi = fix(psi);
            psi = (0:pi/psi:pi)';
        end
        if islogical(lmax) && ~lmax
            spk = false;
        elseif islogical(lmax) && lmax
            spk     = true;
            if psi(1) == 0
                lmax = fix((length(psi) - 1)/2);
            else
                lmax = fix(length(psi)/2);
            end
            wf  = sin(psi)*pi/lmax/2/2;
        elseif ~islogical(lmax) && isscalar(lmax)
            spk = true
            if lmax > fix(length(psi)/2)
                lmax = fix(size(w,1)/2);
                fprintf('WARNING: LMAX cannot be more than half of the number of sampling points when the \n\t sampling is not on the Gauss-Neumann sampling points. \n Taking LMAX = %g \n',lmax)
            end
            wf  = sin(psi)*pi/lmax/2/2;
        else
            error('LMAX must be logical or an integer. Please verify.')
        end
    else
        error('PSI should be a vector or the character string ''Gauss''. Please verify the input.')
    end
elseif nargin == 5 && islogical(spk)
    if ischar(psi) && strcmp(psi,'Gauss')
        [tmp,wf]    = grule(lmax+1);
        wf          = flipud(wf(:))/2;
        psi         = acos(flipud(tmp(:)));
        if isscalar(lmax)
            lmax = fix(lmax);
        else
            error('LMAX must be logical or an integer. Please verify.')
        end
    elseif ~isempty(psi)
        if isvector(psi)
            psi = psi(:);
        elseif isscalar(psi)
            psi = fix(psi);
            psi = (0:pi/psi:pi)';
        end
        if isscalar(lmax)
            if lmax > fix(length(psi)/2)
                lmax = fix(size(w,1)/2);
                fprintf('WARNING: LMAX cannot be more than half of the number of sampling points when the \n\t sampling is not on the Gauss-Neumann sampling points. \n Taking LMAX = %g \n',lmax)
            else
                lmax = fix(lmax);
            end
            wf  = sin(psi)*pi/lmax/2/2;
        else
            error('LMAX must be logical or an integer. Please verify.')
        end
    else
        error('PSI should be a vector or the character string ''Gauss''. Please verify the input.')
    end
else 
    error('Please verify your inputs.')
end

w = 1./sqrt(1+(psi./psi0).^(2*k));

if spk
    P = legpol_rad((0:lmax)',psi);
    w = w.*wf;
    W = P'*w;

    varargout{1} = W;
    varargout{2} = w;
    varargout{3} = psi;
else
    varargout{1} = w;
    varargout{2} = psi;
end
    

