function [mn,del,crms] = monthlyres(ts,mfx)

% MONTHLYRES computes the monthly residuals of a time-series after removing
% a month-wise mean.
%
% [mn,del,crms] = monthlyres(ts)
% [mn,del,crms] = monthlyres(ts,mfx)
%
% INPUT
% ts    -   Time series [t x catchments]. First two columns must be [year
%           month]
% mfx   -   Month-fix switch. Fills the missing months with NaN values.
%
% OUTPUT
% mn    -   Monthly mean
% del   -   Monthly residuals
% crms  -   Catchment-wise RMS value
%--------------------------------------------------------------------------
% USES monthfix
%--------------------------------------------------------------------------

% 16 March 2009, Stuttgart
% Rev. 7 August 2012, Stuttgart 
%   - [monthfix optional]
%--------------------------------------------------------------------------

if nargin == 1
    mfx = false;
elseif nargin == 2
    mfx = logical(mfx);
end

mnth = unique(ts(:,2));
m = length(mnth);

del = ts;
mn = [mnth zeros(m,size(ts,2)-2)];
for i = 1:m
    tmp = del(del(:,2)==mnth(i),3:end);
    N = isnan(tmp);
    n = sum(N);
    tmp(N) = 0;
    tmpmn = (sum(tmp)./(size(tmp,1)-n));
    tmp = tmp - ones(size(tmp,1),1)*tmpmn;
    tmp(N) = NaN;
    del(del(:,2)==mnth(i),3:end) = tmp;
    mn(i,2:end) = tmpmn;
end

n = sum(isnan(del(:,3:end)));
N = find(isnan(del));
del(N) = 0;
crms = sqrt(sum(del(:,3:end).^2)./(size(del,1)-n))';
del(N) = NaN;

if mfx
    del = monthfix(del);
end
