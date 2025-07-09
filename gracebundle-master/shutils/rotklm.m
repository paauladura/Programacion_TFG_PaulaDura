function [rklm,D] = rotklm(klm,lmax,clat,lam)

% ROTKLM rotates the Spherical harmonic spectrum of a field from Earth-centric
% frame to topo-centric frame.
%
% NOTE: This is a function file written specifically for the FilterBundle, and
% so it cannot be used in general for all rotation purposes.
%
% [rklm,D] = rotklm(klm,lmax,clat,lam)
%
% INPUT
% klm   -   complex spherical harmonic coefficients in CS, SC, or
%           look-up-table formats.
% lmax  -   Maximum degree of spherical harmonic expansion.
% clat  -   Co-latitude of the position in degrees.
% lam   -   Longitude of the position in degrees.
%
% OUTPUT
% rklm  -   Rotated complex spherical harmonic coefficients.
% D     -   Rotation matrix or the Wigner D-matrix.
%--------------------------------------------------------------------------
% USES cs2sc sc2cs gcoef2sc cssc2clm dlmk
%--------------------------------------------------------------------------
% See also rotspkcov wignerD
%--------------------------------------------------------------------------
%
%

% 31 May 2011, Stuttgart. Balaji Devaraju
%--------------------------------------------------------------------------

if nargin < 2
    error('Insufficient input arguments')
elseif nargin == 2
    clat = pi/6;
    lam  = pi/6;
elseif nargin == 3
    lam = 0;
end

[m,n] = size(klm);

if m==n && m==lmax+1
    rl  = cs2sc(real(klm),0);
    im  = cs2sc(imag(klm),0);
    klm = rl + 1i*im;
elseif m==((lmax+1)*(lmax+2)/2) && n==4
    rl  = gcoef2sc([klm(:,1:2) real(klm(:,3:4))]);
    im  = gcoef2sc([klm(:,1:2) imag(klm(:,3:4))]);
    klm = rl + 1i*im;
elseif m~=lmax+1 || n~=(2*lmax + 1)
    error('Input SH coefficients not in the required format')
end

c = sum(1:lmax+1);
s = sum(1:lmax);

if lmax <= 100
    D    = struct('cc',zeros(c),'cs',zeros(c,s),'sc',zeros(s,c),'ss',zeros(s));
    D.cc = mat2cell(D.cc,1:lmax+1,1:lmax+1);
    D.cs = mat2cell(D.cs,1:lmax+1,1:lmax);
    D.sc = mat2cell(D.sc,1:lmax,1:lmax+1);
    D.ss = mat2cell(D.ss,1:lmax,1:lmax);
else
    D = 'Warning: (lmax+1)^4 elements cannot be stored in the memory';
end

rklm = klm;

for k = 1:lmax
    tmp = diag(exp(1i*(-k:k)*lam)) * dlmk(k,-clat);
    idx = lmax+1-k:lmax+1+k;
    sc  = klm(k+1,idx);
    
    rklm(k+1,idx) = sc*tmp;
    
    if lmax <= 100
        D.cc{k+1,k+1} = tmp(k+1:end,k+1:end);
        D.ss{k,k}     = rot90(tmp(1:k,1:k),2);
        D.sc{k,k+1}   = flipud(tmp(1:k,k+1:end));
        D.cs{k+1,k}   = fliplr(tmp(k+1:end,1:k));
    end
end

if lmax<=100
    D = [cell2mat(D.cc) cell2mat(D.cs); cell2mat(D.sc) cell2mat(D.ss)];
end
