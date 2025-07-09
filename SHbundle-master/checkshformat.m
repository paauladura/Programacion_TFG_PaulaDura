function [frmt, lmax] = checkshformat(sh)

% CHECKSHFORMAT checks the format of the spherical harmonic
% coefficients.
% 
% [frmt, lmax] = checkshformat(sh)
% 
% INPUT
% sh        - Spherical harmonics in one of the valid formats
%
% OUTPUT
% frmt  - A string denoting one of the valid formats.
% lmax  - Maximum degree of spherical harmonic expansion
%
% The valid formats are
%   'sc'        - /S|C\ spherical harmonic coefficients arranged in a
%                 matrix of size [L+1 x 2L+1]
%                    0   0   0  C00  0   0   0
%                    0   0  S11 C10 C11  0   0
%                    0  S22 S21 C20 C21 C22  0
%                   S33 S32 S31 C30 C31 C32 C33
%   'cs'        - C\S spherical harmonic coefficients arranged in a
%                 matrix of size [L+1 x L+1]
%                   C00 S11 S21 S31
%                   C10 C11 S22 S32
%                   C20 C21 C22 S33
%                   C30 C31 C32 C33
%   'clm_deg'   - Look up table format: degree-leading ordering
%                   l m Clm Slm
%                   -----------
%                   0 0 C00 S00
%                   1 0 C10 S10
%                   1 1 C11 S11
%                   2 0 C20 S20
%                   2 1 C21 S21
%                   2 2 C22 S22
%   'clm_ord'   - Look up table format: order-leading ordering
%                   l m Clm Slm
%                   -----------
%                   0 0 C00 S00
%                   1 0 C10 S10
%                   2 0 C20 S20
%                   1 1 C11 S11
%                   2 1 C21 S21
%                   2 2 C22 S22
%   'klm'       - Look up table format: Degree-wise ordering
%                   l   m   klm
%                   -----------
%                   0   0   C00
%                   1   -1  S11
%                   1   0   C10
%                   1   1   C11
%                   2   -2  S22
%                   2   -1  S21
%                   2   0   C20
%                   2   1   C21
%                   2   2   C22
%
%-----------------------------------------------------------------------



% 
%  Project: Earth GRAvity Field from Space (EGRAFS)
%  Copyright 2015-16 Balaji Devaraju (BD) 
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
%   BD      2016:11:27  Initial version
%----------------------------------------------------------------------

[r, c] = size(sh);

if r == c
    frmt = 'cs';
    lmax = r-1;
elseif c == (2*(r-1) + 1)
    frmt = 'sc';
    lmax = r-1;
elseif c == 4
    l = sh(:,1);
    m = sh(:,2);
    if isequal(sh(:,1), sort(l))
        frmt = 'clm_deg';
        lmax = max(l);
    elseif isequal(sh(:,2), sort(m))
        frmt = 'clm_ord';
        lmax = max(l);
    elseif any(m < 0)
        frmt = 'klm';
        lmax = max(l);
    else
        error('Unknown look-up-table format')
    end
else
    error(sprintf('Invalid format: [%g x %g]', [r, c]))
end



% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
