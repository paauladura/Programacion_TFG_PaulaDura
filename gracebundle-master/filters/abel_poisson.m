function Bl = abel_poisson(h, lmax)

% ABEL_POISSON computes the Abel-Poisson kernel for a given factor h and
% a maximum degree of spherical harmonic expansion
%
% Bl = abel_poisson()
% Bl = abel_poisson(h, lmax)
%
% INPUT
% h     - A factor between 0 and 1 that determines the width of the
%         kernel. A factor of 1 generates the Shannon kernel/box-car
%         function. [default: 0.5]
% lmax  - Maximum degree of spherical harmonic expansion [default: 60]
%
% OUTPUT
% Bl    - Spectral weights of the Abel-Poisson kernel
%
%-----------------------------------------------------------------------



% 
%  Project: GRACE Bundle
%  http://gracebundle.tuxfamily.org
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
%   BD      2017:09:04  Initial version
%----------------------------------------------------------------------

narginchk(0,2)


if nargin==0
    h = 0.5;
    lmax = 60;
elseif nargin==1
    if h > 1
        error('Value of h is greater than 1. Its value should be 0 < h <= 1')
    end
    lmax = 60;
else
    if h > 1
        error('Value of h is greater than 1. Expecting a value of 0 < h <= 1')
    end
    if ~isscalar(lmax)
        error('LMAX is not a scalar. Expecting an scalar integer.')
    end
end


Bl = h.^(0:lmax)';

    




% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
