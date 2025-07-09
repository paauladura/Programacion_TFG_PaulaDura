function W = spk2spt(w,lmax,psi)

% SPK2SPT converts isotropic spectral filters to their spatial
% co-efficients using Legengre polynomials.
%
% This numerical integration technique follows the following convention for
% Legendre polynomial transform
%   W(psi) = sum_0^L ( (2*l+1) * w_l * P_l(cos(psi)) )
%   w_l = 1/2 int_0^pi ( W(psi) P_l(cos(psi)) sin(psi) dpsi )
%
% W = spk2spt(w,lmax)
% W = spk2spt(w,lmax,psi)
%
% INPUT
% w     - spectral weights
% lmax  - maximum degree of spherical harmonic expansion
% psi   - spherical distance in radians
%
% OUTPUT
% W 	- spatial weights [Sph. dist.(radians) Weights] (if 'psi' not given),
%         else only the weights will be given
%------------------------------------------------------------------------------

% USES	legpol

% Balaji Devaraju. Stuttgart 5 May 2009
%------------------------------------------------------------------------------

if nargin == 2
    psi = (0:0.25:180)'*pi/180;
elseif nargin == 3
    psi = psi(:);
else
    error('Insufficient input arguments')
end
n = (0:lmax)';
w = (2*n + 1).*w(:);
P = legpol(0:lmax,psi);

W = P*w;
if nargin == 2
    W = [psi W];
end
