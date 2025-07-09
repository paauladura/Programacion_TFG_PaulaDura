function W = spkcosine(lcut,k,lstart,lmax)

% SPKCOSINE computes the Legendre coefficients of a spectral cosine filter.
% This is the spectral analog of the von Hann filter defined on the sphere,
% and like the von Hann filter this is also a homogeneous isotropic filter.
%
% W = spkcosine
% W = spkcosine(lcut,k,lstart,lmax)
%
% INPUT
% lcut  -   Cut-off degree of the filter [def. lcut = 60]
% k     -   Order of the cosine function. The von Hann is a second order
%           cosine filter, i.e., it is the square of the cosine function.
%           This must be a positive integer, and for very large numbers the
%           coefficients degenerate. [def. k = 2]
% lstart-   The starting point for the cosine taper [def. lstart = 15]
% lmax  -   Maximum degree of Legendre expansion. [def. lmax = lcut]
%
% OUTPUT
% W     -   Legendre coefficients of the cosine filter of order 'k'
%
%--------------------------------------------------------------------------

% Authors:
% 21 April 2011, Stuttgart. Balaji Devaraju.

% Revision history
% BD    01/04/2014  -- Rearranged inputs and brushed up code and comments
%--------------------------------------------------------------------------

if nargin == 0
    lcut    = 60;
    k       = 2;
    lstart  = 15;
    lmax    = 60;
elseif nargin == 1
    k       = 2;
    lstart  = 0;
    lmax    = lcut;
elseif nargin == 2
    lstart  = 0;
    lmax    = lcut;
elseif nargin == 3
    lmax    = lcut;
elseif nargin > 4
    error('Too many input arguments')
end

W           = (0:lmax)';
W(W<lstart) = 1;
W(W>lcut)   = 0;
ind         = ((W ~= 0) & (W >= lstart));
W(ind)      = (cos((W(ind) - lstart)'/(lcut-lstart) * pi/2)).^k;
W(1)        = 1;
