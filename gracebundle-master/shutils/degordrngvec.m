function [klmnew,varargout] = degordrngvec(klm)

% DEGORDRNGVEC orders a given vector of spherical harmonic degrees in
% ascending order. This function is particularly useful if one wants to
% convert the given geopotential coefficient set from Colombo ordering to
% degree ordering.
%
% klmnew               = degordrngvec(klm)
% [klmnew,ind,mind]    = degordrngvec(klm)
%
% INPUT
% klm - Degree vector [deg]/Geopotential coefficient set given in a four 
%        column vector [deg ord Clm Slm] and ordered in Colombo ordering 
%        scheme.
%
% OUTPUT
% klmnew    - Degree ordered vector
% ind       - Index of degree ordering
% mind      - vector containing zeros and ones. Ones indicate non-zero
%              orders. This can be used for selecting Slm coefficients of
%              orders other than zero.
%--------------------------------------------------------------------------

% Created on: 29 March 2008, Stuttgart.
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

[klmnew,ind] = sort(klm(:,1));
klmnew = klm(ind,:);

if nargout==2
    varargout{1} = ind;
elseif nargout == 3
    varargout{1} = ind;
    varargout{2} = (klmnew(:,2)~=0);
end
