function [da,c] = dealias(a, qa)

% DEALIAS removes the aliased tidal residual signal from the residuals of
% the monthly GRACE coefficients/mass estimates.
%
% da        = dealias(a)
% [da, c]   = dealias(a, qa)
%
% INPUT
% a     - Time series with the aliased tidal residuals in it. The data 
%         must be in the following format
%			[year month start_day end_day data]
% qa    - Standard deviation associated with the data
%
% OUTPUT
% da    - Dealiased time series signal.
% c     - Fourier coefficients of the aliased frequencies
%--------------------------------------------------------------------------
% USES grctimetag
%--------------------------------------------------------------------------

% Balaji Devaraju, 9 September 2011, Stuttgart

% Updates
% 2014-01-08 BD Included primary and secondary aliasing signals

if nargin == 0
    error('Insufficient input arguments')
elseif nargin == 1
    qa = [];
elseif nargin == 2
    if ~isempty(qa) && ~isequal(a(:,1:2), qa(:,1:2)) || ~isequal(size(a,2),size(qa,2))
        error(['Mismatch between the time-series of observations', ... '
                'and their standard deviations'])
    end
end

[t,st,et] = grctimetag(a(:,1),a(:,3),a(:,4));

F = [ones(length(st),1) ... % a0
    cos(2*pi*st/160.8) sin(2*pi*st/160.8) ...   % S2 Primary tide aliasing
    cos(2*pi*st/321.7) sin(2*pi*st/321.7) ...   % S1
    cos(2*pi*st/169.73) sin(2*pi*st/169.73) ... % P1
    cos(2*pi*st/2443.8) sin(2*pi*st/2443.8) ... % K1
    cos(2*pi*st/1618.5) sin(2*pi*st/1618.5) ... % K2
    cos(2*pi*st/125.15) sin(2*pi*st/125.15) ...   % M2 Secondary tide aliasing
    cos(2*pi*st/89.03) sin(2*pi*st/89.03) ...   % N2
    cos(2*pi*st/130.8) sin(2*pi*st/130.8) ... % O1
    cos(2*pi*st/90.9) sin(2*pi*st/90.9) ... % Q1
    ];

% N  = F'*F;              % Computing normal matrix
% c  = inv(N)*(F'*a(:,3:end)); % Computing the coefficients

c = zeros(size(F,2),size(a,2)-4);

if isempty(qa)
    c = F\a(:,5:end);
else
    for k = 5:size(a,2)
        P       = 1./(qa(:,k).^2);
        N       = bsxfun(@times,F,P);
        y       = N'*a(:,k);
        N       = N'*F;
        c(:,k)  = N\y;
    end
end

% Removing the tidal alias signal
da          = a;
da(:,5:end) = a(:,5:end) - F(:,2:end)*c(2:end,:); 
