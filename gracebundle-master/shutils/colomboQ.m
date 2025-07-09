function Q = colomboQ(Q,L)

% COLOMBOQ orders a varaiance-covariance matrix, ordered in degree 
% ordering, in Colombo/order-wise ordering.
%
% Q = colomboQ(Q,lmax)
%
% INPUT
% Q     - Full covaraince matrix in degree-ordered format.
%           The filename also could be provided where the string also 
%           includes the path of the file.
% lmax  - Maximum degree of spherical harmonic expansion.
%
% OUTPUT
% Q - Colombo ordered varaince-covariance matrix.
%
% See also colombo, cssc2clm , vcm2vec, degordrngQ
%--------------------------------------------------------------------------
% USES cssc2clm
%--------------------------------------------------------------------------

% Author: Balaji Devaraju
% Created on 18 June 2008, Stuttgart
%--------------------------------------------------------------------------

if nargin < 2
    error('Insufficient inputs')
end

if ischar(Q)
    Q = load(Q);
    tmp = fieldnames(Q);
    Q = Q.(tmp{1});
end

if size(Q,1)~=(L+1)^2 || size(Q,1)~=size(Q,2)
    error('Covariance matrix is either not square or does not contain all the elements')
end

% [r,c] = size(Q);

tmp = cssc2clm([(0:L)' ones(L+1,1)],L);
[t,ic] = sort(tmp(:,1));
tmp = tmp(ic,:);
[t,ic] = sort(tmp(:,2));
[t,is] = sort(tmp((tmp(:,2)~=0),2));
is = is+sum(1:L+1);

% tmp = zeros(r,c);

tmp = Q([ic;is],:);
Q = tmp(:,[ic;is]);
