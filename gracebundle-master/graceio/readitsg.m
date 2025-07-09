function gfc = readitsg(data_path, extn)


% READITSG reads the spherical harmonic coefficients of the ITSG GRACE
% solutions provided by the Institute of Theoretical and Satellite 
% Geodesy, TU Graz. The coefficients are given in ICGEM format, and this
% function acts as a wrapper around the ICGEMPARSER function to extract
% all the coefficients. The output is the 10-column cell array format
%
% gfc = readitsg(data_path, extn)
%
% INPUT
% data_path - Path where all the data is located
% extn      - Extension of the individual ICGEM files
%
% OUTPUT
% gfc   - 10-column cell array format
%           [organisation, coefficient_type, release_number, year, ...
%               month, [start_day end_day total_days], flag, ...
%               max_degree, C\S, standard_deviation]
%
%-----------------------------------------------------------------------
% USES: icgemparser
%-----------------------------------------------------------------------


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
%   BD      2017:09:27  Initial version
%----------------------------------------------------------------------

narginchk(2, 2)

if ~ischar(data_path) || ~isequal(exist(data_path, 'dir'), 7)
    error('Invalid data path')
end


if ~ischar(extn)
    error('Extension is not a character array')
end

if (isunix || ismac)
    if ~strcmp('/', data_path(end))
        data_path = [data_path, '/'];
    end
elseif ispc
    if ~strcmp('\', data_path(end))
        data_path = [data_path, '\'];
    end
end

fprintf('Getting the list of files with %s extension\n', extn)
list_ = dir([data_path, '*.', lower(extn)]);

fprintf('Reading ITSG GRACE L2 data ...\n')
gfc = cell(length(list_), 10);
for k = 1:length(list_)
    tmp = icgemparser([data_path, list_(k).name]);
    gfc{k, 1}   = 'ITSG';
    gfc{k, 2}   = 'GSM';

    [s, e]      = regexp(list_(k).name, ['_n', num2str(tmp.max_degree), '_']);
    gfc{k, 3}   = str2num(list_(k).name(s-4:s-1));
    gfc{k, 4}   = str2num(list_(k).name(e+1:e+4));
    gfc{k, 5}   = str2num(list_(k).name(e+6:e+7));
    gfc{k, 7}   = 0;
    gfc{k, 8}   = tmp.max_degree;
    gfc{k, 9}   = tmp.gfc.knm;
    gfc{k, 10}  = tmp.gfc.dknm;
end
y = vertcat(gfc{:,4});
m = vertcat(gfc{:,5});
d = daysinmonth(y, m);

tmp = [dayofyear(y, m, ones(size(y))), dayofyear(y, m, d), d];
gfc(:, 6) = mat2cell(tmp, ones(length(y), 1), 3);






% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
