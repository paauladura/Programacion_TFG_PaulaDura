function [spvar,rlen,ormat] = hisospvar(B,lmax)

% HISOSPVAR computes the spatial variance of the given homogeneous isotropic 
% filter. Further, it also computes the metrics of resultant length and 
% orientation matrix.
%
% INPUT
% B     -   Spectrum of the homogeneous isotropic filter a column vector
% lmax  -   Maximum degree of spherical harmonic expansion. It must be a
%           positive integer.
% 
% OUTPUT
% spvar -   Spatial variance of the filter
% rlen  -   Resultant length
% ormat -   Orientation matrix
%
%------------------------------------------------------------------------------

% Balaji Devaraju. Stuttgart, 19 July 2014

if nargin < 2
    error('Insufficient input arguments')
end

if ~isint(lmax)
    lmax = fix(lmax);
end

[r,c] = size(B);

if ~isequal(r,lmax+1) || ~isequal(c,1)
    error('Incorrect input format for the spectral coefficients')
end

[psi,wt]    = grule(181);
psi         = acos(flipud(psi(:)));
wt          = flipud(wt(:));

b = spk2spt(B,lmax,psi);

bsqr = b.^2 .* wt;
totE = sum(bsqr);

ormat       = zeros(3);
ormat(1,1)  = sum(sin(psi).^2 .* bsqr)/totE;
ormat(2,2)  = ormat(1,1);
ormat(3,3)  = sum(cos(psi).^2 .* bsqr)/totE;

spvar = asin(sqrt(ormat(1,1)));
rlen  = acos(sum((cos(psi) .* bsqr)/totE));
