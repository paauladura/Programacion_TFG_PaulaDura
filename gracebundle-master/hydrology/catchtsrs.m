function [ts,cid] = catchtsrs(cts,mfix)

% CATCHTSRS generates the time-series of storage changes in catchments,
% which have been aggregated by the MASSESTMTR function. 
%
% [ts,cid] = catchtsrs(cts)
%
% INPUT
% cts   -   Cell of mass estimates as output by the MASSESTMTR function
% mfix  -   '1' if the missing months must be filled up with NaN, else '0'
% 
% OUTPUT
% ts    -   Time-series of all the catchments included in the catchment
%           index file
% cid   -   Catchment ID
%--------------------------------------------------------------------------
% USES GRACE/monthfix
%
% See also massestmtr
%--------------------------------------------------------------------------

% Balaji Devaraju, 26 March 2012, Stuttgart.
%--------------------------------------------------------------------------


if nargin == 1
    mfix = false;
else 
    mfix = logical(mfix);
end

m = size(cts,1);
n = size(cts{1,6},1);

cid = cts{1,6}(:,1)';
ts  = cell2mat(cts(:,1:3));
ts  = [ts zeros(m,n)];

for i = 1:m
    ts(i,6:end) = cts{i,6}(:,2)';
end

if mfix, ts = monthfix(ts); end

ts = [zeros(1,5) cid; ts];
