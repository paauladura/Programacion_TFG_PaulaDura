function [d,sd,ed] = grctimetag(y,sd,ed)

% GRCTIMETAG generates a time tag for the GRACE dataset whose units are
% days instead of the months. Every GRACE monthly solution comes with the
% information of the days in the year that were used to generate the
% solution, and this information is used for generating the time tag. The
% time tag is the mid-point between the start day and the end day, and
% these time tags are counted continuously.
%
% d         = grctimetag(y,sd,ed)
% [d,sd,ed] = grctimetag(y,sd,ed)
%
% INPUT
% y     -   vector of years of all the monthly solutions
% sd    -   starting day of the monthly solution
% ed    -   end day of the monthly solution
% 
% OUTPUT
% d     -   vector of time tags
% sd    -   starting day of the monthly solution in terms of continuous 
%           days from first day of the first monthly solution.
% ed    -   end day of the monthly solution in the same way as 'sd'.
%--------------------------------------------------------------------------

% Balaji Devaraju. Stuttgart, 26 July 2012.

if nargin < 3
    error('Insufficient input arguments')
end

if length(y(:))~=length(sd(:)) || length(sd(:))~=length(ed(:))
    error('The size of the inputs do not match')
end

uy  = unique(y);
diy = daysinyear(uy-1);

diy(uy==uy(1)) = 0;
diy = cumsum(diy);

ey = y;
ey(sd>ed) = ey(sd>ed) + 1;

ds = y;
de = ey;
for k = 1:length(uy)
    ds(y==uy(k)) = diy(k);
    de(ey==uy(k)) = diy(k);
end
sd = ds+sd;
ed = de+ed;
d  = (sd + ed)/2;
