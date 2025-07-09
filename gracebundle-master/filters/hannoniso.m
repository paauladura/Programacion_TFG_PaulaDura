function B = hannoniso(psi0,mc,psi1,lmax,wndw)

% HANNONISO computes the spectral co-efficient values of the Non-isotropic
% filter developed by SC Han et al, GJI (2005) 163,18--25.
%
% B = hannoniso
% B = hannoniso(psi0)
% B = hannoniso(psi0,wndw)
% B = hannoniso(psi0,mc,psi1,lmax,wndw)
%
% INPUT
% psi0  -   fundamental radius [radians]
% mc    -   order threshold
% psi1  -   secondary radius [radians] (optional) [default: psi1 = psi0*2]
% lmax  -   maximum degree of the SH development
% wndw  -   Type of smoothing window to be applied. (optional)
%           'gauss' - Gaussian smoothing window [default]
%           'hann'  - Hann smoothing window
%           'pell'  - Pellinen smoothing window
%
% OUTPUT
% B     -   spectral co-efficients of the filter [l m Clm Slm] given in
%           Colombo ordering.
%--------------------------------------------------------------------------

% USES
% 	FilterBundle/gaussfltr
% 	             vonhann
% 	             pellinen
%   SHbundle/cssc2clm

% Created on 19 September 2007
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

if nargin == 0
    psi0    = pi/36;
    mc      = 36;
    psi1    = psi0*2;
    lmax    = 3*mc;
    wndw    = 'gauss';
elseif nargin == 1
    mc      = fix(pi/psi0);
    psi1    = psi0*2;
    lmax    = 3*mc;
    wndw    = 'gauss';
elseif nargin == 2
    if ~ischar(mc)
        wndw = 'gauss';
    elseif ischar(mc)
        wndw    = mc;
        mc      = pi/psi0;
    end
    psi1 = psi0*2;
    lmax = 3*mc;
elseif nargin == 3
    if ~ischar(psi1)
        wndw = 'gauss';
    elseif ischar(psi1)
        wndw = psi1;
        psi1 = psi0*2;
    end
    if ~isint(mc)
        mc = fix(mc);
    end
    lmax = 3*mc;
elseif nargin == 4
    if ~isint(mc)
        mc = fix(mc);
    end
    if ~ischar(lmax)
        wndw = 'gauss';
        if ~isint(lmax)
            lmax = fix(lmax);
        end
    elseif ischar(lmax)
        wndw = lmax;
        lmax = 3*mc;
    end
elseif nargin == 5
    if ~isint(mc)
        mc = fix(mc);
    end
    if ~isint(lmax)
        lmax = fix(lmax);
    end
end

B = cssc2clm(zeros(lmax+1),lmax);

for m = 1:lmax+1
    temp    = [(m-1:lmax)',ones(lmax-m+2,1).*(m-1)];
    psi2    = ((psi1-psi0)/mc * (m-1)) + psi0;
    ind     = (B(:,2)==(m-1));
    if strcmp('gauss',wndw)
        Btmp        = gaussfltr(psi2,lmax);
        B(ind,3)    = Btmp(m:end,1);
    elseif strcmp('hann',wndw)
        Btmp        = vonhann(psi2*2,lmax);
        B(ind,3)    = Btmp(m:end);
    elseif strcmp('pell',wndw)
        B(ind,3) = pellinen((m-1:lmax)',(psi2*2));
    else
        error('Unknown smoothing window specified')
    end
end
B(:,4) = B(:,3);
