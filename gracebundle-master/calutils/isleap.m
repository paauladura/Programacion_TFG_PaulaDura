function idx = isleap(y)

% ISLEAP determines if a given year or a vector of years is a leap year or not.
% 
% idx = isleap(y)
%
% INPUT
% y     - vector of years
%
% OUTPUT
% idx   - Same size as y, with logical indices: 'true' if leap year, else 
%         'false'.
%------------------------------------------------------------------------------

% Balaji Devaraju, Stuttgart, 1 October 2012
%------------------------------------------------------------------------------


idx = false(size(y));

idx(mod(y,400)==0) = true;
idx(mod(y,100)~=0 & mod(y,4)==0) = true;
