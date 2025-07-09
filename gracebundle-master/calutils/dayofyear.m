function n = dayofyear(y,m,d)

% DAYOFYEAR calculates the day of the year for a given date
%
% n = dayofyear(y,m,d)
%
% INPUT
% y - year
% m - month
% d - day
% All inputs must be of the same size
%
% OUTPUT
% n - day of the year
%--------------------------------------------------------------------------

% Created on: 3 September 2009, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

if nargin<3 || nargin>3
    error('Insufficient or too many input arguments')
elseif ~isequal(length(y),length(m)) || ~isequal(length(m),length(d))
    error('Check the length of the date vectors')
end

leap        = cumsum([0;31;29;31;30;31;30;31;31;30;31;30]);
nleap       = cumsum([0;31;28;31;30;31;30;31;31;30;31;30]);

n 		    = zeros(size(y));
icent 	    = (mod(y,400)==0);
icent 	    = icent + (mod(y,100)~=0 & mod(y,4)==0);
icent 	    = logical(icent);

n(icent)	= leap(m(icent)) + d(icent);
n(~icent)	= nleap(m(~icent)) + d(~icent);

n((m(~icent) == 2) & (d(~icent) > 28)) = NaN;
n((m(icent)  == 2) & (d(icent)  > 29)) = NaN;
n((d < 1) | (d > 31)) = NaN;
n((m < 1) | (m > 12)) = NaN;
