function beta = pellinen(psi,l)

% PELLINEN(psi,l) returns Pellinen's smoothing beta_n factors.
%
% beta = pellinen(psi)
% beta = pellinen(psi,l)
% 
% INPUT
% psi - radius of smoothing cap [radians]		vector
% l   - degree                                  vector
%
% OUTPUT
% beta - Pellinen beta factors.
%       If both L and PSI are vectors, L indicates the columns, and PSI the rows
%       If one (or both) of them is scalar, the output vector follows the shape 
%       of the input.  
%-------------------------------------------------------------------------------

% Nico Sneeuw                        Munich                           17/08/94
%
% uses LEGPOL_RAD
% rev.1 27/09/96 NS: 
%    - proper handling of l=0 and of psi=0
%    - a bit brushing up
%       15/02/2013
%    - input psi changed to radians
%    - legpol --> legpol_rad


% diagnostics 
if nargin < 2, 
    l = (0:60); 
end

if min(size(l)) > 1 || min(size(psi)) > 1, 
    error('L and PSI should be vectors (or scalars).'); 
end

if any(psi < 0) || any(psi > pi), 
    error('PSI must be within [0;pi] degrees.'); 
end

if any(l< 0),
    error('The minimum degree should be non-negative.'), 
end

% preliminaries
[lrow,lcol] = size(l);
[prow,pcol] = size(psi);
l           = l(:)';					% row vector
psi         = psi(:);					% column "
lmax        = max(l);
lmin        = min(l);


% Create matrix of unnormalized Legendre functions, and difference the
% proper columns. Create columns of ones (no zeros!) in case lmin=0.
if lmin == 0
   pmat = legpol([0 0:lmax+1],psi);
else
   pmat = legpol(lmin-1:lmax+1,psi);
end
lind = l - lmin + 2;				% index into columns

% Calculate the beta's. Special care for psi=0.
beta        = pmat(:,lind-1) - pmat(:,lind+1);
m           = (psi==0);					% mask for psi=0
beta(m,:)   = 1;
if sum(m) < length(psi)				% i.e. non-zero psi exists
   beta(~m,:) = beta(~m,:) .* ( (1./(1-cos(psi(~m)))) * (1./(2*l+1)) );
end

% reshape beta if either l or psi is scalar
if numel(l) == 1,   beta = reshape(beta,prow,pcol); end
if numel(psi) == 1, beta = reshape(beta,lrow,lcol); end
