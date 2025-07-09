function [W,D] = rotspkcov(W,lmax,clat,lam,intyp)

%
% ROTSPKCOV rotates spectral covariance/smoothing functions on the sphere
% using the Euler rotations (alpha,beta,gamma). The rotations are
% accomplished by the use of Wigner-D matrices. These rotations are valid
% only for normalized covariance/smoothing function matrices, which have
% been normalized using the geodetic conventions.
%
% W    = rotspkcov(W,lmax,clat,lam,intyp)
%[W,D] = rotspkcov(W,lmax,clat,lam,intyp)
%
% INPUT
% W     -   Spectral coefficient matrix of covariance/smoothing functions
%           in degree-leading format. The matrix must be square and
%           complex.
% lmax  -   Maximum degree of spherical harmonic expansion.
% lam   -   Longitude of the location of the covariance function.
% clat  -   Co-latitude of the location of the covariance function.
% intyp -   Input structure of the W matrix. 
%           'full'  - Fully populated W matrix
%           'iso'   - Isotropic functions given as [l Wl]
%           'cs'    - Diagonal W-matrix stored in CS-format.
%
% OUTPUT
% W     -   Rotated covariance/smoothing function matrix
% D     -   Wigner-D matrix used for rotating W-matrix
%--------------------------------------------------------------------------
%
% USES dlmk cs2sc
%--------------------------------------------------------------------------
%
% See also rotklm
%--------------------------------------------------------------------------
%
%

if nargin < 5
    error('Insufficient input arguments')
end

if mod(lmax,floor(lmax))~=0
    error('lmax must be an integere')
end

if strcmp(intyp,'full')
    if size(W,1)~=size(W,2) || size(W,1)~=(lmax+1)^2
        error('Check the dimensions of input spectral covariance/smoothing function matrix')
    end
    
    c = sum(1:lmax+1);
    s = sum(1:lmax);
    
    D    = struct('cc',zeros(c),'cs',zeros(c,s),'sc',zeros(s,c),'ss',zeros(s));
    D.cc = mat2cell(D.cc,1:lmax+1,1:lmax+1);
    D.cs = mat2cell(D.cs,1:lmax+1,1:lmax);
    D.sc = mat2cell(D.sc,1:lmax,1:lmax+1);
    D.ss = mat2cell(D.ss,1:lmax,1:lmax);
    
    D.cc{1} = 1;
    for k = 1:lmax
        tmp = diag(exp(1i*(-k:k)'*lam))* dlmk(k,-clat); % Computation of rotation matrix
                
        D.cc{k+1,k+1} = tmp(k+1:end,k+1:end);
        D.ss{k,k}     = rot90(tmp(1:k,1:k),2);
        D.sc{k,k+1}   = flipud(tmp(1:k,k+1:end));
        D.cs{k+1,k}   = fliplr(tmp(k+1:end,1:k));
    end
    D = [cell2mat(D.cc) cell2mat(D.cs); cell2mat(D.sc) cell2mat(D.ss)];
    W = W*D; % Rotation of coefficients
elseif strcmp(intyp,'iso')
    if size(W,1) ~= lmax+1 || size(W,2) ~= 2
        error('Check the dimensions of the input spectrum of isotropic function. W must be of the form [l Wl]')
    end
    
    c = sum(1:lmax+1);
    s = sum(1:lmax);
    
    D    = struct('cc',zeros(c),'cs',zeros(c,s),'sc',zeros(s,c),'ss',zeros(s));
    D.cc = mat2cell(D.cc,1:lmax+1,1:lmax+1);
    D.cs = mat2cell(D.cs,1:lmax+1,1:lmax);
    D.sc = mat2cell(D.sc,1:lmax,1:lmax+1);
    D.ss = mat2cell(D.ss,1:lmax,1:lmax);
    
    temp = D;
    
    D.cc{1} = 1;
    temp.cc{1} = W(1,2)*D.cc{1};
    for k = 1:lmax
        tmp = diag(exp(1i*(-k:k)'*lam))* dlmk(k,-clat); % Computation of rotation matrix
        
        D.cc{k+1,k+1} = tmp(k+1:end,k+1:end);
        D.ss{k,k}     = rot90(tmp(1:k,1:k),2);
        D.sc{k,k+1}   = flipud(tmp(1:k,k+1:end));
        D.cs{k+1,k}   = fliplr(tmp(k+1:end,1:k));
        
        tmp = W(k+1,2)*tmp; % Rotation of coefficients
        temp.cc{k+1,k+1} = tmp(k+1:end,k+1:end);
        temp.ss{k,k}     = rot90(tmp(1:k,1:k),2);
        temp.sc{k,k+1}   = flipud(tmp(1:k,k+1:end));
        temp.cs{k+1,k}   = fliplr(tmp(k+1:end,1:k));
    end
    D = [cell2mat(D.cc) cell2mat(D.cs); cell2mat(D.sc) cell2mat(D.ss)];
    W = [cell2mat(temp.cc) cell2mat(temp.cs); cell2mat(temp.sc) cell2mat(temp.ss)];
    
elseif strcmp(intyp,'cs')
    if size(W,1)~=size(W,2) || size(W,1)~=lmax+1
        error('Check the dimensions of input spectral covariance/smoothing function matrix')
    end

    c = sum(1:lmax+1);
    s = sum(1:lmax);
    
    D    = struct('cc',zeros(c),'cs',zeros(c,s),'sc',zeros(s,c),'ss',zeros(s));
    D.cc = mat2cell(D.cc,1:lmax+1,1:lmax+1);
    D.cs = mat2cell(D.cs,1:lmax+1,1:lmax);
    D.sc = mat2cell(D.sc,1:lmax,1:lmax+1);
    D.ss = mat2cell(D.ss,1:lmax,1:lmax);
    
    sc = cs2sc(W);
     W = D;
    
    D.cc{1} = 1;
    W.cc{1} = sc(1,lmax+1);
    
    for k = 1:lmax
        tmp = diag(exp(1i*(-k:k)'*lam))* dlmk(k,-clat); % Computation of rotation matrix
        
        D.cc{k+1,k+1} = tmp(k+1:end,k+1:end);
        D.ss{k,k}     = rot90(tmp(1:k,1:k),2);
        D.sc{k,k+1}   = flipud(tmp(1:k,k+1:end));
        D.cs{k+1,k}   = fliplr(tmp(k+1:end,1:k));
        
        tmp = diag(sc(k+1,lmax+1-k:lmax+1+k)) * tmp; % rotation of coefficients
        W.cc{k+1,k+1} = tmp(k+1:end,k+1:end);
        W.ss{k,k}     = rot90(tmp(1:k,1:k),2);
        W.sc{k,k+1}   = flipud(tmp(1:k,k+1:end));
        W.cs{k+1,k}   = fliplr(tmp(k+1:end,1:k));
    end
    D = [cell2mat(D.cc) cell2mat(D.cs); cell2mat(D.sc) cell2mat(D.ss)];
    W = [cell2mat(W.cc) cell2mat(W.cs); cell2mat(W.sc) cell2mat(W.ss)];
end
