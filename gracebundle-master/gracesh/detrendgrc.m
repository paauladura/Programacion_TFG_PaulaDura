function [dts,tnd,A,Qx] = detrendgrc(ts,qts,n)

% DETREND estimates and removes a trend from a given monthly time-series
% dataset.
%
% [dts,tnd,A,Qx] = detrend(ts,qts,n)
%
% INPUT
% ts    -   Time-series matrix with the rows denoting time and the columns
%           denoting the observations [t x m]. The format of the matrix must 
%           be as follows: 
%               [year month start_day end_day y1 y2 y3 ...].
% qts   -   Standard deviations of the observations in the same format as
%           the time-series matrix.
%               [year month start_day end_day qy1 qy2 qy3 ...].
% n     -   Fitting a trend is in general fitting a polynomial, and hence,
%           the desired polynomial degree.
%
% OUTPUT
% dts   -   Detrended time-series
% tnd   -   Trend values for each of the m observation points.
% A     -   Design matrix used for estimating trends.
% Qx    -   Cell array containing the covariance matrix of the trend
%           estimates.
%
%--------------------------------------------------------------------------
% USES GRACEBundle/ monthfix
%                   grctimetag
%--------------------------------------------------------------------------

% Balaji Devaraju, 7 September 2011, Stuttgart.

% Checking input arguments
if nargin == 1
    qts = [];
    n   = 1;
elseif nargin == 2
    if ~isempty(qts) && ~isequal(ts(:,1:2), qts(:,1:2))
        error('Temporal mismatch between the time-series of observations and their standard deviations')
    end
    n = 1;
elseif nargin == 3
    if ~isempty(qts) && ~isequal(ts(:,1:2), qts(:,1:2))
        error('Temporal mismatch between the time-series of observations and their standard deviations')
    end
    if ~isint(n)
        fprintf('WARNING: Polynomial degree was not an integer. Using only the integer part. \n')
        n = fix(n);
    end
end

% Initializing values and selecting integral number of years for
% trend estimation.
ts      = monthfix(ts);
t       = (1:size(ts,1))';
t       = t(~isnan(ts(:,3)));
ts      = ts(~isnan(ts(:,3)),:);
T       = find(mod(t,12)==0,1,'last');
fprintf('T = %g \n',T)
[t,st]  = grctimetag(ts(:,1),ts(:,3),ts(:,4));
t       = st;
A       = ones(length(t),n+1);
A(:,2)  = t;
if n > 1
    for k = 3:n+1
        A(:,k) = t.^(k-1);
    end
end

m   = size(ts,2)-4;
tnd = zeros(n+1,m);
Qx  = cell(m,1);

% Computing trend with standard deviations of GRACE coefficients
if ~isempty(qts)
    for k = 1:m
        P           = 1./qts(1:T,k+4).^2;
        N           = bsxfun(@times,A(1:T,:),P);
        y           = N'*ts(1:T,k+4);
        N           = N'*A(1:T,:);
        Qx{k}       = N\eye(size(N));
        x           = Qx{k}*y;
        tnd(:,k)    = x; %.*(abs(x) > 2*(sqrt(diag(Qx{k}))));
    end
% Computing trend without standard deviations of GRACE coefficients
elseif isempty(qts)
    N   = A(1:T,:)'*A(1:T,:);
    y   = A(1:T,:)'*ts(1:T,5:end);
    Q   = N\eye(size(N));
    tnd = Q * y;
    s   = sum((ts(1:T,5:end) - A(1:T,:)*tnd).^2)/(T-n-1);
    for k = 1:m
        Qx{k} = Q*s(k);
    end
end
% Removing the estimated trend from the time-series
dts             = ts;
dts(:,5:end)    = ts(:,5:end) - A(:,2:end)*tnd(2:end,:);
