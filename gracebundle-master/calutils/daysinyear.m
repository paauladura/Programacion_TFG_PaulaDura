function d = daysinyear(y)

% DAYSINYEAR computes the number of days in a year
%
% d = daysinyear(y)
%
% INPUT
% y - Year
%
% OUTPUT
% d - number of days in the year 'y'
%--------------------------------------------------------------------------

% Balaji Devaraju. Stuttgart, 26 July 2012.

y = y(:);

d = zeros(size(y));

for k = 1:size(y,1)
    if mod(y(k),100) == 0
        if mod(y(k),400) == 0, d(k) = 366; else d(k) = 365; end
    else
        if mod(y(k),4) == 0, d(k) = 366; else d(k) = 365; end
    end
end