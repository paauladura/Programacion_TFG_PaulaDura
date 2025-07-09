function fxdts = monthfix(ts)

% MONTHFIX fixes gaps in a time-series with monthly values. The time-series
% has to be in ascending temporal order. This program should mainly be used
% as a tool to find the data gaps and fill them with NaN values.
% Interpolation will not be done within this function. The utility of this
% program is for plotting time-series showing the data gaps properly.
%
% INPUT
% ts    -   Times series with [time x data]. The first two columns must be
%           [year month] and the month column must be integers from 1 to 12.
% 
% OUTPUT
% fxdts -   Time-series that is fixed for temporal data gaps. The data gaps
%           are filled with NaN values.
%--------------------------------------------------------------------------

% Balaji Devaraju, 8 August 2011, Stuttgart.

% BD    2013-03-22  Completely rewritten without for loops.

if ~isempty(find((ts(:,2) > 12) | (ts(:,2) < 1)))
	error('Month values must be a positive integer between 1 and 12')
end

yrs = unique(ts(:,1));
syr = yrs(1);

indx = (ts(:,1)-syr)*12 + ts(:,2);

y = ones(12,1) * yrs';
y = y(:);
m = (1:12)' * ones(1,length(yrs));
m = m(:);

fxdts 			= ones(length(y),size(ts,2))*NaN;
fxdts(:,1) 		= y;
fxdts(:,2) 		= m;
fxdts(indx,:) 	= ts;
fxdts 			= fxdts(indx(1):indx(end),:);
