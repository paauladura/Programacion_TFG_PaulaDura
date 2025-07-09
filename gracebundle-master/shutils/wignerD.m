function D = wignerD(l,a,b,g,opt)

% WIGNERD computes the Wigner-D matrices for a given degree of
% spherical harmonic expansion. These matrices are useful in rotating
% spherical harmonic functions as well as spherical harmonic coefficients.
% The Wigner-D matrices computed here are normalized in the geodetic sense
% as given in Inclination functions, Sneeuw(1991). Reports from TU Delft.
% This function can calculate both real as well as complex Wigner-D
% matrices.
%
% D = wignerD
% D = wignerD(l)
% D = wignerD(l,alpha)
% D = wignerD(l,alpha,beta)
% D = wignerD(l,alpha,beta,gamma)
% D = wignerD(l,alpha,beta,gamma,opt)
%
% INPUT
% l     -   Maximum degree of spherical harmonic expansion.
% a,b,g -   Euler angles a(lpha), b(eta), g(amma) [radians]
%           0 <= a <= 2*pi,0 <= g <= 2*pi and 0 <= b <= pi
% opt   -   Option for computing real ('r') or complex ('c') Wigner-D
%           matrices
%
% OUTPUT
% D     -   Wigner-D matrix for the particular degree
%--------------------------------------------------------------------------
%
% USES dlmk, cpx2realmat
%
% See also wigner2
%--------------------------------------------------------------------------
%
%

% Created on: 14 July 2010
% Author: Balaji Devaraju (BD)
%
% Revision history
% 2014-03-08    BD  Changed 'lmax' to 'l'
%                   Fine-tuned code for efficient computation
%--------------------------------------------------------------------------

if nargin == 0
    l = 3;
    a = pi/2;
    b = a;
    g = a;
    opt = 'c';
elseif nargin == 1
    a = pi/2;
    b = a;
    g = a;
    opt = 'c';
elseif nargin == 2
    b = 0;
    g = 0;
    opt = 'c';
    if a >= 2*pi
        a = a - floor(a/(2*pi)) * 2*pi;
    elseif a < 0 && a >= -pi
        a = pi-a;
    elseif a < -pi
        error('Check the value of ''a'' it must be 0 < a < 360')
    end
elseif nargin == 3
    g = 0;
    opt = 'c';
elseif nargin == 4
    opt = 'c';
elseif nargin == 5
    if isempty(a), a = pi/2; end
    if isempty(b), b = pi/2; end
    if isempty(g), g = pi/2; end
end

if ~isint(l)
    l = floor(l);
end

m = -l:l;
D = (exp(-1i*m'*g) * exp(-1i*m*a)).*dlmk(l,b);
if strcmp(opt,'r')
    C = cpx2realmat(l);
    D = C*D*C';
    D = real(D);
end
D = D';
