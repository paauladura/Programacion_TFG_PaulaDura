function klm_slr = cellreplacec20(klm, slr)

%
% CELLREPLACEC20 replaces C20 value of GRACE using SLR values provided
% via Technical Note 7
%
% klm_slr = cellreplacec20(klm, slr)
%
% INPUT
% klm   - Accepts 10-column cell output given by READGRC or the time-
%         series matrix output by KLM2TSRS
% slr   - SLR C20 values as output by RDSLRC20
%
% OUTPUT
% klm_slr - Cell-array/Matrix with the replaced C20
%
%-----------------------------------------------------------------------
% See also rdslrc20
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

narginchk(2, 2)

if iscell(klm) && (size(klm, 2) < 9)
    error('Expecting a 10-column cell array')
end

klm_slr  = klm;
t = [vertcat(klm_slr{:, 4}), vertcat(klm_slr{:, 6})];

for k = 1:length(t)
    idx = ((slr(:,3)==t(k,1)) & (abs(slr(:,4)-t(k,2))<=10));
    if any(idx)
        klm_slr{k,9}(3,1)   = slr(idx,5);
        klm_slr{k,10}(3,1)  = slr(idx,7);
    else
        klm_slr{k,9}(3,1)   = NaN;
        klm_slr{k,10}(3,1)  = NaN;
    end
end


% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
