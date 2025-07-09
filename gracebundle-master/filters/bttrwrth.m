function bl = bttrwrth(lc,order,lmax)

% BTTRWRTH computes the spectral filter coefficients of a Butterworth
% filter. Butterworth filter is an isotropic filter.
%
% bl = bttrwrth(lc,order,lmax)
% 
% INPUT
% lc    -   Degree after which the expansion must be filtered.
% order -   Power to which the ratio has to be raised
% lmax  -   Maximum degree of spherical harmonic expansion
% 
% OUTPUT
% bl    -   coefficients of the spectral Butterworth filter.
%-----------------------------------------------------------------------
% See also spbttrwrth
%-----------------------------------------------------------------------
% uses
% 	cssc2clm
%-----------------------------------------------------------------------


% Created on: 1 March 2008, Stuttgart
% Author: Balaji Devaraju
%-----------------------------------------------------------------------

if nargin == 3
    if lmax < 0 || lc < 0 || order < 0
        error('One of the inputs is negative. Inputs must be positive')
    end
    if lmax < lc
        error('LMAX is greater than the cut-off degree')
    end
elseif nargin == 2
    if order < 0 || lc < 0
        error('One of the inputs is negative. Inputs must be positive')
    end
    lmax    = 2 * lc;
elseif nargin == 1
    if lc < 0 
        error('One of the inputs is negative. Inputs must be positive')
    end
    lmax    = 2 * lc;
    order 	= 2;
elseif nargin == 0
    lmax    = 120
    lc      = 40
    order   = 2;
end

bl = (0:lmax)';
bl = 1./sqrt(1+(bl./lc).^(2*order));

