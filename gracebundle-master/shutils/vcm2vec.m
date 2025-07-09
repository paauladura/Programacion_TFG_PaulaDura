function dvlist = vcm2vec(Qxx, lmax)

% VCM2VEC converts the diagonal elements of the variance-covariance matrix
% of a given set of spherical harmonic co-efficients into the
% Colombo-ordering, i.e., [l m(sorted) Clm Slm]
%
% dvlist = vcm2vec(Qxx,lmax)
%
% INPUT 
% Qxx 	- Variance-Covariance Matrix
% lmax 	- Maximum degree of the Spherical Harmonic co-efficients
%
% OUTPUT
% dvlist - [l m(sorted) Sigma_Clm^2 Sigma_Slm^2]
%--------------------------------------------------------------------------

% Created on 30 May 2007
% Author: Balaji Devaraju, Stuttgart
%--------------------------------------------------------------------------

lmcount = sum(1:lmax+1);

l = zeros(lmcount,1);
m = zeros(lmcount,1);
cnt = 1;
for i = 0:lmax
    m(cnt:i+cnt) = (0:i)';
    l(cnt:i+cnt) = ones(i+1,1).*i;
    cnt = cnt+i+1;
end

[ms,indm] = sort(m);

var = diag(Qxx);
varc = var(1:lmcount);
vars = zeros(lmcount,1);
vars(m~=0) = var(lmcount+1:end);

dvlist = [l(indm) ms varc(indm) vars(indm)];
