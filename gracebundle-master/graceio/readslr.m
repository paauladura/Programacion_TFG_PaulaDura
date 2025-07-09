function gfc = readslr(data_path, extn)


% READSLR reads the spherical harmonic coefficients of the SLR
% solutions provided by Astronomical Institute, Bern. The coefficients 
% are given in ICGEM format, and this function acts as a wrapper around 
% the ICGEMPARSER function to extract all the coefficients. The output 
% is the 10-column cell array format.
%
% gfc = readslr(data_path, extn)
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
%   BD      2019:07:10  Initial version
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
    endif
end

fprintf('Getting the list of files with %s extension\n', extn)
list_ = dir([data_path, '*.', lower(extn)]);

fprintf('Reading hl-SST L2 data ...\n')
gfc = cell(length(list_), 10);
for k = 1:length(list_)
    tmp = icgemparser([data_path, list_(k).name]);
    
    gfc{k, 1}   = 'AIUB';
    gfc{k, 2}   = 'GSM';
    gfc{k, 3}   = 0;
    gfc{k, 4}   = str2num(tmp.time_period{1}(1:4));
    gfc{k, 5}   = str2num(tmp.time_period{1}(5:6));
    start_day   = dayofyear(gfc{k,4}, gfc{k,5}, str2num(tmp.time_period{1}(7:8)));
    end_day     = dayofyear(str2num(tmp.time_period{1}(10:13)), ...
                            str2num(tmp.time_period{1}(14:15)), ...
                            str2num(tmp.time_period{1}(16:17)));
    gfc{k, 6}   = [start_day, end_day, (end_day - start_day + 1)];
    gfc{k, 7}   = 0;
    
    
    gfc{k, 8}   = tmp.max_degree;
    gfc{k, 9}   = tmp.gfc.knm;
    gfc{k, 10}  = tmp.gfc.dknm;
end



% vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
