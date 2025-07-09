function cs = klmtsrs2sc(klm,ordrng,flg)

% KLMTSRS2SC converts rows of spherical harmonic coefficients from KLMTSRSMAT
% to CS/SC formats. This is complementary function for the KLMTSRSMAT function.
%
% cs = klmtsrs2sc(klm_ts,ordrng)
%
% I/P
% klm_ts    - spherical harmonic coefficients from KLMTSRSMAT. The format
%             of this matrix must be [year month start_day end_day data]
% ordrng    - Ordering of the spherical harmonic coefficients as given by 
%             KLMTSRSMAT function.
% flg       - Flag for indicating SC ('0') and CS ('1') formats.
%
% O/P
% cs        - spherical harmonic coefficients in SC/CS format, with each row
%             in the klm matrix stored as a four column cell array.
%               { year month [start_day end_day] cs }
%------------------------------------------------------------------------------
% USES clm2sc, sc2cs
%
% See also klmtsrsmat
%------------------------------------------------------------------------------

% Balaji Devaraju, Stuttgart, 10 January 2010
%
% Revision
% 2014-01-14 BD Added flag for SC and CS format choice
%------------------------------------------------------------------------------

if nargin < 2
    error('Insufficient input arguments')
elseif nargin == 2
    flg = 0;
end
flg = logical(flg);

lmax    = max(ordrng(:,1));
cff     = (lmax+1)^2;

[m,n] = size(klm);

if m ~= (cff + 4) && n == (cff + 4)
    klm = klm';
    n   = size(klm,2);
elseif n ~= (cff + 4)
    error('Complete set of spherical harmonic coefficients not provided')
end

cs  = [mat2cell(klm(1:4,:)',ones(n,1),[1,1,2]), cell(n,1)];
klm = klm(5:end,:);

c   = sum(1:lmax+1);
for k = 1:n
    tmp = [ordrng(1:c,:) klm(1:c,k) zeros(c,1)];
    tmp(tmp(:,2)~=0,4) = klm(c+1:end,k)';
    tmp = clm2sc(tmp);
    if flg
        tmp = sc2cs(tmp);
    end
    cs{k,4} = tmp;
end
