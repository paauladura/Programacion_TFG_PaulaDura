function [s, e] = doymonth(y, m)

% DAYOFYEARMONTH computes the day of year of the beginning and end of a
% given month with year.
%
% [s, e] = doymonth(y, m)
%
% INPUT
% y - year  [n x 1]
% m - month [n x 1]
%
% OUTPUT
% s - day of year of the staring of the month   [n x 1]
% e - day of year of the end of the month       [n x 1]
%
%-----------------------------------------------------------------------

% Created on: 3 August 2016, Hannover
% Author: Balaji Devaraju
%-----------------------------------------------------------------------

narginchk(2, 2)

% Check the dimensions of the inputs
[r1, c1] = size(y);
[r2, c2] = size(m);

if (r1 == r2) && (c1 == c2)
    if r1 == 1
        y   = y';
        m   = m';
    end
elseif (r1==c2) && (c1 == r2)
    m = m';
else
    error(sprintf(['Dimension mismatch between year and day of the year variables.\n', ...
                    'Year and day of the year variables must be of the same size.']))
end

% Is it a leap year?
icent 	  = (mod(y,400)==0);
icent 	  = icent + (mod(y,100)~=0 & mod(y,4)==0);
icent 	  = logical(icent);

% Start day of the month
s 		  = zeros(size(y));
s_leap    = cumsum([1;31;29;31;30;31;30;31;31;30;31;30]);
s_nleap   = cumsum([1;31;28;31;30;31;30;31;31;30;31;30]);
s(icent)  = s_leap(m(icent));
s(~icent) = s_nleap(m(~icent));
s((m < 1) || (m > 12)) = NaN;

% end day of the month
e         = s;
e_leap    = cumsum([31;29;31;30;31;30;31;31;30;31;30;31]);
e_nleap   = cumsum([31;28;31;30;31;30;31;31;30;31;30;31]);
e(icent)  = e_leap(m(icent));
e(~icent) = e_nleap(m(~icent));
e((m < 1) || (m > 12)) = NaN;
