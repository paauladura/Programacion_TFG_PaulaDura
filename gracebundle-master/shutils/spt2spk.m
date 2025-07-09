function W = spt2spk(w,psi,lmax)

% SPT2SPK converts an isotropic function defined on a sphere into its
% spectrum using Legendre polynomials
%
% This numerical integration technique follows the following convention for
% Legendre polynomial transform
%   w(psi) = sum_0^L ( (2*l+1) * W_l * P_l(cos(psi)) )
%   W_l    = 1/2 int_0^pi ( w(psi) P_l(cos(psi)) sin(psi) dpsi )
%
% W = spt2spk(w,psi)
% W = spt2spk(w,psi,lmax)
%
% INPUT
% w   	- spatial function values [mx1] or [mxn]
% psi 	- Sampling points in radians. If the spatial function was sampled on a
%         Gauss-Neumann grid then specify 'Gauss' instead of the values of the
%         sampling points.
% lmax 	- Degree of maximum expansion of the Legendre polynomials. Not required
%         if the second input argument is 'Gauss'.
%
% OUTPUT
% W   	- spectral coefficients
%--------------------------------------------------------------------------

% USES	legpol_rad

% Balaji Devaraju. Stuttgart 5 May 2009
%--------------------------------------------------------------------------

if nargin < 2
    error('Insufficient input arguments')
elseif nargin == 2
    if isvector(w)
        w = w(:);
    end
    if ischar(psi) && strcmp(psi,'Gauss')
        [tmp,wf]    = grule(size(w,1));
        wf          = flipud(wf(:));
        psi         = acos(flipud(tmp(:)));
        lmax        = size(w,1)-1;
        G           = true;
    elseif isequal(length(psi),size(w,1))
        if psi(1) == 0
            lmax = fix((size(w,1) - 1)/2);
        else
            lmax = fix(size(w,1)/2);
        end
        G = false;
    else
        error('Check input variable PSI')
    end
elseif nargin == 3
    if size(w,1)==1 && size(w,2)>1
        w = w(:);
    end
    if ~isint(lmax)
        lmax = fix(lmax);
    end
    if lmax > fix(size(w,1)/2)
        lmax = fix(size(w,1)/2);
        fprintf('WARNING: LMAX cannot be more than half of the number of sampling points when the \n\t sampling is not on the Gauss-Neumann sampling points. \n Taking LMAX = %g \n',lmax)
    end

    G = false;
end

psi = psi(:);
L 	= (0:lmax)';
P   = legpol_rad(L,psi);

if G
    wf  = wf/2;
else
    wf  = sin(psi)*(pi/lmax/2)/2;
end
w = bsxfun(@times,w,wf);
W = P'*w;
