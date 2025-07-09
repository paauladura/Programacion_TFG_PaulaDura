function [m, d] = doy2cal(y, dy)

% DOY2CAL outputs the calendar date given the year and day of the year.
%
% INPUT
% y     - year              [1 x 1]
% dy    - day of the year   [1 x 1]
% All inputs must be of the same size.
%
% OUTPUT
% m, d  - Month and day of the month. [n x 1]
%
%-----------------------------------------------------------------------

% Created on: 3 August 2016, Hannover
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

% Check the number of input arguments
narginchk(2, 2)

% Check the dimensions of the inputs
[r1, c1] = size(y);
[r2, c2] = size(dy);

if (r1 == r2) && (c1 == c2)
    if r1 == 1
        y   = y';
        dy  = dy';
    end
elseif (r1==c2) && (c1 == r2)
    dy = dy';
else
    error(sprintf(['Dimension mismatch between year and day of the year variables.\n', ...
                    'Year and day of the year variables must be of the same size.']))
end

[m, d] = arrayfun(@day2cal_main, y, dy, 'UniformOutput', true);

end


function [m, d] = day2cal_main(y, dy)

leap        = cumsum([0;31;29;31;30;31;30;31;31;30;31;30]);
nleap       = cumsum([0;31;28;31;30;31;30;31;31;30;31;30]);

n 		    = zeros(size(y));
icent 	    = (mod(y,400)==0);
icent 	    = icent + (mod(y,100)~=0 & mod(y,4)==0);
icent 	    = logical(icent);

if icent
    m = find(ismember(leap, dy));
    if isempty(m)
        dff = dy - leap;
        d = min(dff(dff > 0));
        m = find(dff == d);
    else
        m = m - 1;
        d = dy - leap(m);
    end
else
    m = find(ismember(nleap, dy));
    if isempty(m)
        dff = dy - nleap;
        d = min(dff(dff > 0));
        m = find(dff == d);
    else
        m = m - 1;
        d = dy - nleap(m);
    end
end

end
