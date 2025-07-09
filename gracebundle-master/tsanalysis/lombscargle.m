function [f,P,sgf,sigma] = lombscargle(x,y,ofac,fmax)

% LOMBSCARGLE computes the power spectral density of a given
% one-dimensional dataset that is unevenly sampled using the Lomb-Scargle
% algorithm. This algorithm is explained thoroughly in the references given
% below. Most of the variable names have been retained as given in Press et
% al. (2007), but this implementation is a vectorized direct implementation
% of the formulas given there.
% References:
% Press et al. 2007. Numerical Recipes §13.8
% Horne and Baliunas, 1986. The Astrophysical Journal, 302 (2):757-763.
%
% INPUT
% x     -   sampling points.
% y     -   observed data. This can be a column vector or a matrix. If it
%           is a matrix then the operations are performed column-wise.
% ofac  -   scaling factor that decides the spacing for the frequency
%           search.
% fmax  -   highest frequency (not angular) that needs to be searched for.
%           If no input is given then Nyquist frequency is chosen as the
%           'fmax'. This must have same units as the sampling points.
%
% OUTPUT
% f     -   vector of frequencies whose power was computed.
% P     -   estimated power of the frequencies in 'f'. P has the same 
%           number of columns as the observed dataset 'y'.
% sgf   -   significance of the power of the each of the estimated
%           frequencies 'f'.
% sigma -   computed satndard deviation of the dataset 'y'.
%
%--------------------------------------------------------------------------

% Balaji Devaraju. Stuttgart, 26 July 2012.

if length(x) < 2
    error('Insufficient observations to carry out periodogram analysis')
end

if size(y,1) ~= length(x)
    error('Sizes of x and y do not match')
end

[n,m] = size(y);
xdiff = max(x) - min(x);
xave  = (max(x) + min(x))/2;
fc    = n/(2*xdiff);

if nargin == 2
    ofac = 4;
    fmax = fc;
elseif nargin == 3
    fmax = fc;
end

hifac = fmax/fc;
np = floor(ofac*hifac*n/2);

f   = (1:np)'/(xdiff*ofac); % frequencies whose power will be estimated
arg = 2*pi*f*(x-xave)';

yave  = mean(y);            % Mean of y
sigma = std(y);             % Standard deviation
dy    = y - ones(n,1)*yave; % Residual

Sh = sin(arg)*dy;       % trigonometric recursions
Ch = cos(arg)*dy;
S2 = sum(sin(2*arg),2);
C2 = sum(cos(2*arg),2);

tau  = atan2(S2,C2); % tau here is actually 2*2*pi*f*tau

ctau = cos(tau/2) * ones(1,m); % Numerator
stau = sin(tau/2) * ones(1,m);
Cnum = Ch.*ctau + Sh.*stau;
Snum = Sh.*ctau - Ch.*stau;

CS2  = C2.*cos(tau) + S2.*sin(tau); % Denominator
Cden = (n + CS2)/2;
Sden = (n - CS2)/2;

P   = Cnum.^2./(Cden*2*(sigma.^2)) + Snum.^2./(Sden*2*(sigma.^2)); % Periodogram
sgf = 1 - (1 - exp(-P)).^(2*np/ofac); % significance
