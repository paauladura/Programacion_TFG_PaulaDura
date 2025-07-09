function Q = degordrngQ(Q,L)

% DEGORDRNGQ orders a varaiance-covariance matrix, ordered in order-leading 
% ordering, in degree-leading ordering.
%
% Q = degordrngQ(Q,lmax)
%
% INPUT
% Q     - Full covaraince matrix in order-leading format.
%           The filename also could be provided where the string also 
%           includes the path of the file.
% lmax  - Maximum degree of spherical harmonic expansion.
%
% OUTPUT
% Q - Degree ordered varaince-covariance matrix.
%
% See also degordrngvec, cssc2clm, vcm2vec, colomboQ
%--------------------------------------------------------------------------
% USES cssc2clm
%--------------------------------------------------------------------------

% Author: Balaji Devaraju
% Created on 28 March 2009, Stuttgart
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

tmp = cssc2clm([(0:L)' ones(L+1,1)],L);
[t,ic] = sort(tmp(:,1));
[t,is] = sort(tmp((tmp(:,2)~=0),1));
is = is+sum(1:L+1);

tmp = Q([ic;is],:);
Q = tmp(:,[ic;is]);
