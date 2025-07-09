function klm_slr = tsreplacec20(klm, idx, slr, typ)

%
% TSREPLACEC20 replaces C20 value of GRACE using SLR values provided
% via Technical Note 7
%
% klm_slr = tsreplacec20(klm, idx, slr)
%
% INPUT
% klm - The time-series matrix output by KLM2TSRS
%           [year month start_day end_day klm]
% idx - Column index of the C20 coefficient
% slr - SLR C20 values as output by RDSLRC20
% typ - Whether the time-series contains 'estimate', or 'std'
%       (standard deviation)
%
% OUTPUT
% klm_slr - Matrix with the replaced C20
%
%-----------------------------------------------------------------------
% See also rdslrc20, cellreplacec20
%-----------------------------------------------------------------------
%

% Balaji Devaraju
% 
%  Project: gracebundle
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
%   BD      2017:09:18  Initial version
%----------------------------------------------------------------------

narginchk(3, 4)

if (nargin==4) && ~ischar(typ)
    error('Type of GRACE input must be a character array')
end

if ~isscalar(idx)
    error('Index for C20 coefficient is not a scalar value')
end

klm_slr  = klm;
t = klm_slr(:, [1, 3, 4]);

switch lower(typ)
case { 'estimate' }
    for k = 1:length(t)
        indx = ((slr(:,3)==t(k,1)) & (abs(slr(:,4)-t(k,2))<=10));
        if any(idx)
            klm_slr(k, idx) = slr(indx,5);
        else
            klm_slr(k, idx) = NaN;
        end
    end
case { 'std' }
    for k = 1:length(t)
        indx = ((slr(:,3)==t(k,1)) & (abs(slr(:,4)-t(k,2))<=10));
        if any(idx)
            klm_slr(k, idx) = slr(indx,7);
        else
            klm_slr(k, idx) = NaN;
        end
    end
end


% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

