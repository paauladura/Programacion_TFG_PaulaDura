function ymdays = daysinmonth(year,month)

% DAYSINMONTH calculates the number of days for a given year and 
% month.
%
% ymdays = daysinmonth(year,month)
% 
% INPUT
% year  - Can be scalar or vector.
% month - Can be a scalar or vector.
%
% If both month and year are vectors they have to be of the same size.
%
% OUTPUT
% ymdays - A vector of days for the given months.
%--------------------------------------------------------------------------

% Created on: 10 Mar 2008
% Author: Balaji Devaraju, Stuttgart
%--------------------------------------------------------------------------

if nargin<2
    error('Insufficient inputs')
end

if any(month<1) || any(month>12)
    error('Invalid inputs for months')
end

if sum(size(year)) == 2
    year = ones(size(month)).*year;
elseif sum(size(month)) == 2
    month = ones(size(year)).*month;
elseif size(year) ~= size(month)
    error('Year has to be either a scalar or has to match the size of the month vector')
end

n = size([year month],1);
ymdays = zeros(n,1);

jan2jul = (month<=7);
aug2dec = (month>7);
oddmnth = (mod(month,2)==1);
evnmnth = (mod(month,2)==0);

ymdays(jan2jul.*oddmnth == 1) = 31;
ymdays(aug2dec.*evnmnth == 1) = 31;
ymdays(jan2jul.*evnmnth == 1) = 30;
ymdays(aug2dec.*oddmnth == 1) = 30;

feb = (month==2);
nrmlyear = (mod(year,4)~=0);
nrmlcen = (mod(year,100)==0 & mod(year,400)~=0);
leapcen = (mod(year,400)==0);
leapyear = (mod(year,100)~=0 & mod(year,4)==0);

ymdays(feb.*leapcen == 1) = 29;
ymdays(feb.*leapyear == 1) = 29;
ymdays(feb.*nrmlyear == 1) = 28;
ymdays(feb.*nrmlcen == 1) = 28;
