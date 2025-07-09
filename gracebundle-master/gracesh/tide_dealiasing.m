function [da, klmdot] = tide_dealiasing(a, t_mean, qa)

% TIDE_DEALIASING performs dealiasing of aliased tidal frequencies
%
% a      - data [year month start_day end_day data]
% t_mean - Temporal origin for computing the trend values [1 x 1] [days]
% qa     - Standard deviation of the data  (optional)
%           [year month start_day end_day std_data]
%
% da - dealiased dataset
% klmdot - trend
%------------------------------------------------------------------------




%
%  Project: GRACE Bundle
%  Copyright Balaji Devaraju (BD)
%  devaraju at ife dot uni-hannover dot de
%
%  License: GNU GPLv3 or later
%  You should have received a copy of the GNU General Public License
%  along with EGRAFS;  If not, see <http://www.gnu.org/licenses/>.
%
%  Authors: Balaji Devaraju
%
%  Version control
%  Auhtor   YYYY:MM:DD  Comment
% 	BD 		2017:09:18 	Initial version
%----------------------------------------------------------------------


% Getting continuous time tags
[~, st, et] = grctimetag(a(:,1), a(:,3), a(:,4));

% Referencing to time origin
t = st - t_mean;

% Design matrix
F = [   ones(length(t),1)                       ... % a0
        t                                       ... % trend
        cos(2*pi*t/182.5)  sin(2*pi*t/182.5)    ... % Semiannual
        cos(2*pi*t/365)    sin(2*pi*t/365)      ... % Annual
        cos(2*pi*t/160.8)  sin(2*pi*t/160.8)    ... % S2 Primary tide aliases
        cos(2*pi*t/321.7)  sin(2*pi*t/321.7)    ... % S1
        cos(2*pi*t/169.73) sin(2*pi*t/169.73)   ... % P1
        cos(2*pi*t/2443.8) sin(2*pi*t/2443.8)   ... % K1
        cos(2*pi*t/1618.5) sin(2*pi*t/1618.5)   ... % K2
        cos(2*pi*t/125.15) sin(2*pi*t/125.15)   ... % M2 Secondary tide aliases
        cos(2*pi*t/89.03)  sin(2*pi*t/89.03)    ... % N2
        cos(2*pi*t/130.8)  sin(2*pi*t/130.8)    ... % O1
        cos(2*pi*t/90.9)   sin(2*pi*t/90.9)     ... % Q1
    ];

% Initializing the coefficients
c = zeros(size(F,2), size(a,2)-4);
if nargin < 3
    c = F\a(:, 5:end);
else
    for k = 5:size(a,2)
        P = 1./(qa(:,k).^2);
        N = bsxfun(@times, F, P);
        y = N' * a(:,k);
        N = N' * F;
        c(:, k-4) = N\y;
    end
end

% Removing the tidal alias signal
da = a;
da(:, 5:end) = a(:, 5:end) - F(:, 7:end)*c(7:end, :);

klmdot = c(2,:);
