function B = vonhann(psi0,lmax)

% VONHANN computes the Hann window spectral coefficients for a given
% cap radius. The formulas are taken from S. Becker, 2004, Dipl. Thesis,
% Uni. Bonn.
%
% B = vonhann(psi0,lmax)
% 
% INPUT
% psi0  -   averaging radius in radians [def: pi/36]
% lmax  -   maximum degree [def: 120]
%
% OUTPUT
% B     -   Spectral coefficients of the von Hann smoothing function.
%
%--------------------------------------------------------------------------

% USES 	legpol_rad

% Created on: 30 November 2007, Stuttgart
% Author: Balaji Devaraju

% Revision History:
% BD    02/04/2014 Bug fixing 

%--------------------------------------------------------------------------

if nargin == 0
    lmax    = 120;
    psi0    = pi/36; % 5 degrees
elseif nargin == 1
    lmax    = 120;
elseif nargin == 2
    if ~isint(lmax)
        lmax = fix(lmax);
        fprintf('WARNING: LMAX was not an integer. Using FIX(LMAX). \n')
    end
    if (psi0 < 0) || (psi0 > pi)
        error('Averaging radius PSI0 must be in radians. Please verify.')
    end
end

B 	= [1; zeros(lmax,1)];
psi = (0:psi0/1e4:psi0)';
Pn  = legpol_rad(0:lmax+1,psi);
Pn 	= Pn';

beta 	= (pi^2 -psi0^2)/(pi^2*(1 - cos(psi0)) - 2*psi0^2);
cnst 	= cos((pi/psi0)*psi).*sin(psi).*(psi0/10000);
gma 	= Pn(2:end-1,:) * cnst;
Pldiff 	= Pn(1:(end-2),end) - Pn(3:end,end);
linv 	= (1./((2*(1:lmax)')+1));

B(2:end,1) = beta.*(linv.*Pldiff + gma);
