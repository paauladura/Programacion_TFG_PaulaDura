function B = gaussfltr(psi0,lmax,typ)

% GAUSSFLTR calculates the Gaussian filter weight co-efficients with the given
% maximum degree and the averaging radius. The coefficients are computed by 
% numerical integration.
%
% B = gaussfltr
% B = gaussfltr(psi0,lmax)
%
% I/P 
% psi0  -   averaging radius in radians [def: pi/36]
% lmax  -   maximum degree [def: 120]
% typ   -   method for computing the coefficients pi/lmax is greater than psi0
%           'alias' - Numerically integrate with lmax samples
%           'trunc' - Numerically integrate with 2*pi/psi0 samples [default]
%
% O/P 
% B     -   Spectral coefficients of the Gaussian smoothing function.
%------------------------------------------------------------------------------
% USES  uberall/grule
%------------------------------------------------------------------------------

% Created on: 20 February 2007, Stuttgart
% Author: Balaji Devaraju
%
% Revision:
% 2014-01-27    BD      Complete re-write. No more recursion formulas. Uses
%                       numerical integration of the analytical formula.
% 2014-10-05    BD      Included variable "typ" and distinguished aliased 
%                       spectrum and truncated spectrum, when "lmax" does not 
%                       provide enough sampling points to sample the filter at 
%                       the prescribed radius.
%------------------------------------------------------------------------------

if nargin == 0
    lmax    = 120;
    psi0    = pi/36; % 5 degrees
    L       = lmax;
    typ     = 'trunc';
elseif nargin == 1
    lmax    = round(2*pi/psi0);
    L       = lmax;
    typ     = 'trunc';
elseif nargin == 2
    if ~isint(lmax)
        lmax = fix(lmax);
        fprintf('WARNING: LMAX was not an integer. Using FIX(LMAX). \n')
    end
    if psi0 < 4*pi/lmax
        L = round(4*pi/psi0);
    else
        L = lmax;
    end
    typ = 'trunc';
    if (psi0 < 0) || (psi0 > pi)
        error('Averaging radius PSI0 must be in radians. Please verify.')
    end
elseif nargin == 3
    if ~isint(lmax)
        lmax = fix(lmax);
        fprintf('WARNING: LMAX was not an integer. Using FIX(LMAX). \n')
    end
    if strcmp(typ,'trunc')
        if psi0 < 4*pi/lmax
            L = round(4*pi/psi0);
        else
            L = lmax;
        end
    elseif strcmp(typ,'alias')
        L = lmax;
    else
        error('TYP variable indiscernible. Please verify.')
    end
end

[psi,wf]    = grule(L+1); % Computing the zeros of the Legendre polynomial of degree lmax.
psi         = flipud(acos(psi'));
wf          = flipud(wf');

a       = log(2)/(1 - cos(psi0));
w       = exp(-a*(1-cos(psi))); % unnormalized filter values
w       = w.*wf;

P       = legpol(0:lmax,psi);

B       = P'*w/2;
B       = B./B(1);
B(B<0)  = 0;


if strcmp(typ,'trunc')
    B = B(1:lmax+1);
end
