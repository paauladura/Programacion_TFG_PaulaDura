function [c, s] = blddsgn(pos,lmax,quant,const)

% BLDDSGN builds a design matrix of solid spherical harmonics given the 
% position of the satellite and the maximum degree of spherical harmonic 
% expansion. The solid spherical harmonics are arranged in order-leading format.
% 
% [c s] = blddsgn(pos,lmax,quant,const)
% 
% INPUT
% pos   -   [nx3] vector with longitude and co-latitude all in radians.
%           Height of the satellite must be given from the center of the
%           earth in [m].
%           [longitude, co-latitude, height] 
% lmax  -   [1x1] scalar integer value of the maximum degree of 
%           spherical harmonic expansion.
% quant -   character string which is any of the following: 'none' 
%           (coefficients define the output), 'geoid', 'potential',
%           'dg' or 'gravity' (grav. anomaly), 'tr' (grav. disturbance),
%           or 'trr' (2nd rad. derivative), 'water' (water equivalent 
%           height)
%                                                           - def: 'none'
% const -   Gravitational constant GM and semi-major axis of the reference
%           ellipsoid in [m]
%           [GM ae] (optional)                      - def: GRS80 constants
%
% OUTPUT
% c,s   -   Clm and Slm part of the design matrix provided separately. The
%           spherical harmonics are arranged in order-leading format.
%--------------------------------------------------------------------------
% USES  SHbundle/cssc2clm, eigengrav, upwcon, plm
%       uberall/constants
%--------------------------------------------------------------------------

% Created on: 5 February 2008, Stuttgart
% Author: Balaji Devaraju
%
% Revision history:
%   2014-03-04  BD  Code optimization and included 'const' variable
%                   brushed up help text
%--------------------------------------------------------------------------

%------------------
% Checking inputs 
%------------------
[prow,pcol] = size(pos);
if nargin == 1
    constants
    const   = [GM ae];
    if pcol < 2
        error('Position data does not conform to requirements')
    elseif pcol == 2
        pos = [pos ones(prow,1).*const(2)];
    end
    lmax    =   70;
    quant   =   'none';
elseif nargin == 2
    constants
    const   = [GM ae];
    if pcol < 2 
        error('Position data does not conform to requirements')
    elseif pcol == 2
        pos = [pos ones(prow,1).*const(2)];
    end
    if (mod((lmax*10),10)~=0)
        lmax = floor(lmax);
        display('Warning: Value provided for maximum degree was not an integer. Maximum degree is now equal to the floored value')
    end
    quant   =   'none';
elseif nargin == 3
    constants
    const   = [GM ae];
    if pcol < 2 
        error('Position data does not conform to requirements')
    elseif pcol == 2
        pos = [pos ones(prow,1).*const(2)];
    end
    if (mod((lmax*10),10)~=0)
        lmax = floor(lmax);
        display('Warning: Value provided for maximum degree was not an integer. Maximum degree is now equal to the floored value')
    end
elseif nargin == 4
    if pcol < 2 
        error('Position data does not conform to requirements')
    elseif pcol == 2
        pos = [pos ones(prow,1).*const(2)];
    end
    if (mod((lmax*10),10)~=0)
        lmax = floor(lmax);
        display('Warning: Value provided for maximum degree was not an integer. Maximum degree is now equal to the floored value')
    end
    if numel(const)>2
        error('Variable CONST must be either [2x1] or [1x2] vector')
    end
end

%------------------------------------------------
% Calculating isotropic transfer coefficients
%------------------------------------------------
transf  = eigengrav((0:lmax)', quant, 0, const);
transf  = cssc2clm([(0:lmax)' transf],lmax);
l       = transf(:,1);
m       = transf(:,2);
transf  = transf(:,3);

%---------------------------------------
% Calculating the Legendre functions
%---------------------------------------
P = zeros(size(pos,1),sum(1:lmax+1));
idx2 = cumsum(lmax+1:-1:1);
idx1 = [1,idx2(1:end-1)+1];
for k = 1:lmax+1
    P(:,idx1(k):idx2(k)) = plm((k-1:lmax),k-1,pos(:,2));
end

%-------------------------------------------------------
% Calculating Clm and Slm parts of the design matrix
%-------------------------------------------------------
c = P.*cos(pos(:,1)*m'); % Clm
P = P(:,lmax+2:end);
m = m(lmax+2:end);
s = P.*sin(pos(:,1)*m'); % Slm
P = []; m = [];

if pcol == 3 % Compute upward continuation only if radius was provided
    transf = (const(2)./pos(:,3))*transf';
    pos(:,3) = pos(:,3) - const(2);
    transf = upwcon(l',pos(:,3),const).*transf;
    c = c.*transf;
    transf = transf(:,lmax+2:end);
    s = s.*transf;
else
    c = bsxfun(@times,c,transf');
    transf = transf(lmax+2:end,1)';
    s = bsxfun(@times,s,transf);
end
