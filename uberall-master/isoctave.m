function lgc = isoctave ()

% ISOCTAVE checks if the current environment is octave. If the 
% environment is Octave then it is true else it is false.
%



% 
%  Project: uberall
%  License: GNU GPLv3 or later
%  You should have received a copy of the GNU General Public License
%  along with EGRAFS;  If not, see <http://www.gnu.org/licenses/>.
%  
%  Authors: Balaji Devaraju
% 
%  Version control
%  BD   2016:03:17  Initial version
%  
%-----------------------------------------------------------------------


persistent x;

if isempty(x)
    x = logical(exist('OCTAVE_VERSION', 'builtin'));
end
lgc = x;


% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
