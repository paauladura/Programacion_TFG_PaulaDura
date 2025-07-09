function [b,nrmf] = sptgauss(psi0,psi,nrm)

% SPTGAUSS calculates the Gaussian averaging operator in the spatial domain
% based on the formulas given by Wahr et al 1998. JGR.
%
% b 		= sptgauss(psi0,psi)
% [b,nrmf]  = sptgauss(psi0,psi,nrm)
%
% INPUT
% psi0  -   Averaging radius of the operator [radians] [def: pi/36]
% psi   -   Radii at which the Gaussian averaging operator values are sought.
%           The input must be given as a column vector [radians] [deg: 1deg grid]
% nrm 	- 	Switch for performing normalization of the coefficients. Takes
% 			values '0' or '1' with '1' switching the normalization on. The
% 			default value is '0'
%
% OUTPUT 
% b    -   (Normalized) weights at the respective spatial radii (nx1)[no units]
%           If PSI is not given as an input then b is a (nx2) vector [psi b]
% nrmf 	- 	Normalizing factor
%
% References:
% Jekeli, C. 1981. "Alternative methods to smooth the Earth's gravity field".
% OSU reports 327. Eqn. no. 61.
%
%-------------------------------------------------------------------------------
%


% Created on: 8 August 2007, Stuttgart
% Author: Balaji Devaraju
%-------------------------------------------------------------------------------

if nargin == 0
    psi0    = pi/36;
    psi     = (0:pi/180:pi)';
    nrm     = false;
elseif nargin == 1
    psi     = (0:pi/180:pi)';
    nrm     = false;
elseif nargin == 2
    if (psi0 < 0) || (psi0 > pi)
        error('Averaging radius PSI0 must be in radians. Please verify.')
    end
    if min(size(psi)) == 1
        psi = psi(:);
    else
        error('Variable PSI must be a row/column vector')
    end
	nrm     = false;
elseif nargin == 3
    if (psi0 < 0) || (psi0 > pi)
        error('Averaging radius PSI0 must be in radians. Please verify.')
    end
    if min(size(psi)) == 1
        psi = psi(:);
    else
        error('Variable PSI must be a row/column vector')
    end
    if (nrm == 0) || (nrm == 1)
        nrm = logical(nrm);
    else
        error('Variable NRM must be either ''0'' or ''1''.')
    end
end

a       = log(2)/(1 - cos(psi0));
nrmf    = (1-exp(-2*a))/(2*a);

if nrm
	b = exp(-a*(1-cos(psi)))/nrmf; % normalized coefficients
else
	b = exp(-a*(1-cos(psi))); % unnormalized coefficients
end

if nargin < 2
    b = [psi b];
end
