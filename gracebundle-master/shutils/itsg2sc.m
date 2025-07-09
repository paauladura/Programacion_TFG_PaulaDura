function sc = itsg2sc(klm)

% ITSG2SC converts coefficients arranged in ITSG format to SC-format
%
% INPUT
% klm - Spherical harmonic coefficients arranged in ITSG format
%
% OUTPUT
% sc  - SPherical harmonic coefficients in /S|C\ triangular format
%

narginchk(1, 1)

lmax = max(klm(:,1));
sc = zeros(lmax+1, 2*lmax+1);

for k = 1:length(klm)
    sc(klm(:,1)+1, lmax+1+klm(:,2)) = klm(:,3);
end
