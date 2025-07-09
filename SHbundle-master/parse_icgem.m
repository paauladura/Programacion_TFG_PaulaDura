function [coeffs, maxGOout, minGOout, info] = parse_icgem(filename, varargin)

% PARSE_ICGEM reads files in the ICGEM format.
%
% IN:
%    filename ....... coefficient path/file name
%    'max_lm' ....... reads coefficients up to a degree = order = max_lm
%                     (default: inf, read all coefficients)
%    'min_lm' ....... omit coefficients below that degree/order
%                     (default: 0)
%    'verbose' ...... if true, dumps additional info on the screen
%                     (default: false)
%    'warning' ...... if true, warning about empty coefficients is shown
%                     (default: true)
% OUT:
%    coeffs ......... 'gfc' coefficients from ICGEM file in clm-format
%    maxGOout ....... maximum degree and order read
%    minGOout ....... minimum degree and order read
%    info ........... if available, provides useful information
%                     from the ICGEM file like
%                     - 'earth_gravity_constant',
%                     - 'radius' etc.
%
% EXAMPLE:
%   parse_icgem('./GRACE/ITG-Grace2010s.gfc'); 
%   parse_icgem('./GRACE/ITG-Grace2010s.gfc', 'max_lm', 20); 
%   parse_icgem('./GRACE/ITG-Grace2010s.gfc', 'max_lm', 20, 'min_lm', 5, 'verbose', true); 
%
% SCREEN DUMP:
%    'comment' lines, as well as unknown line identifiers
%
% USES:
%    uberall/GETOPT

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias ROTH (MR), GIS, Uni Stuttgart 
%    Balaji DEVARAJU (BD), IfE, Leibniz Uni Hannover
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2018-07-25: MA, added column-information for GFZ-data (necessary for MATLAB2012?)
%    2015-04-10: BD, brush up help text
%    2014-10-22: MR, put named parameters functionality into a new function
%                    uberall/GETOPT
%    2014-09-25: MR, rename file to 'parse_icgem' to honour the changes in
%                    I/O-parameters and change of functionality (e.g.
%                    named parameters)
%    2014-01-13: MR, strread --> textscan
%    2013-10-18: MR, add params.min_lm to omit lower degrees, workaround for
%                    incorrect line identifier ' gfc': remove leading
%                    whitespace
%    2013-02-14: MR, further help text brush up
%    2012-09-28: MR, a line in the data file containing only a single space 
%                    gave an error
%    2012-06-25: MR, brush up help text,
%                    return maxGO,
%                    disp(sprintf(...)) --> fprintf(...)
%    2012-02-24: MR, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

%% check parameters
if ~exist(filename, 'file')
    error('file "%s" not found...', filename);
end
fprintf('parsing file "%s"...\n', filename);

% define defaults and check optional parameters
defaultparams = {'max_lm',  inf;
                 'min_lm',  0;
                 'verbose', false;
                 'warning', true};
params = getopt(defaultparams, false, varargin);                 

% status output
if params.verbose
    fprintf('  ... trying to read coefficients starting with degree/order "%d"\n', params.min_lm);
    fprintf('  ... trying to read coefficients up to degree/order "%d"\n', params.max_lm);
end

%% Initialize
minGOout = [];
maxGOout = [];
info     = [];

%% First, parse for line identifier and rest of line
lines = {};

fid = fopen(filename, 'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end   
    if isempty(tline), continue, end
    % extract line identifier
    w = 1; % remove leading whitespace
    while ((tline(w) == ' ') || (cast(tline(w), 'uint8') == 9)) && (w < length(tline)) % blank or tab
        w = w + 1;
    end
    i = w;
    while ((tline(i) ~= ' ') && (cast(tline(i), 'uint8') ~= 9)) && (i < length(tline)) 
        i = i + 1;
    end
    j = i;
    while ((tline(j) == ' ') || (cast(tline(j), 'uint8') == 9))
        j = j + 1;
        if length(tline) < j, break; end;
    end;
    
    lines(end + 1, :) = [{tline(w:i-1)}, {tline(j:end)}];
end
fclose(fid);


%% Second, take line identifiers and parse through rest of line accordingly
header = true; % we expect a header first!

err_hand = 1; % Standard: mit sigma C/S einlesen (wird nötig für DDK-Daten)
minGO = params.max_lm;
maxGO = params.min_lm;
tmp_coeffs = cell((maxGO+1)^2, maxGO);
for i = 1:length(lines);
    lineid   = upper(lines{i, 1});
    linecont = lines{i, 2};
    switch lineid
        case 'COMMENT' % comments
            if params.verbose, fprintf('  NOTICE:  %s %s\n', lineid, linecont); end;
        case 'ERRORS'
            err = textscan(linecont, '%s'); 
            if params.verbose, fprintf('  NOTICE:  %s %s\n', lineid, err{1}{1}); end;
            if strcmp(err{1}{1}, 'no')
                err_hand = 0;
            elseif strcmp(err{1}{1}, 'formal')
                err_hand = 1; 
            end
        case 'GFC' % coefficients
            if err_hand == 0
                lin_tmp = textscan(linecont, '%f %f %f %f'); 
                cols = 4;
            elseif err_hand == 1
                lin_tmp = textscan(linecont, '%f %f %f %f %f %f');
                cols = 6;
            end
            if (lin_tmp{1} >= params.min_lm) && (lin_tmp{1} <= params.max_lm) && (lin_tmp{2} <= params.max_lm) % only coeffs minGO <= degree/order <= maxGO
                minGO = min([minGO, lin_tmp{1}]);
                maxGO = max([maxGO, lin_tmp{1}, lin_tmp{2}]);
                % compute correct line for order and degree
                lin = (lin_tmp{1} + 2) * (lin_tmp{1} + 1) / 2 - (lin_tmp{1} - lin_tmp{2});
                tmp_coeffs(lin, 1:cols) = lin_tmp;             
            end
        case 'PRODUCT_TYPE'
            if nargout >= 4
                if params.verbose,  fprintf('  %s %s\n', lineid, linecont); end;
                info.product_type = linecont; 
            else
                if params.verbose, fprintf('  NOTICE:  skipping line %d - no output variable declared: %s %s\n', i, lineid, linecont); end;
            end
        case 'MODELNAME'
            if nargout >= 4
                if params.verbose, fprintf('  %s %s\n', lineid, linecont); end;
                info.modelname = linecont; 
            else
                if params.verbose, fprintf('  NOTICE:  skipping line %d - no output variable declared: %s %s\n', i, lineid, linecont); end;
            end
        case 'RADIUS'
            if nargout >= 4
                if params.verbose, fprintf('  %s %s\n', lineid, linecont); end;
                info.R = cell2mat(textscan(linecont, '%f')); 
            else
                if params.verbose, fprintf('  NOTICE:  skipping line %d - no output variable declared: %s %s\n', i, lineid, linecont); end;
            end
        case 'EARTH_GRAVITY_CONSTANT'
            if nargout >= 4
                if params.verbose, fprintf('  %s %s\n', lineid, linecont); end;
                info.GM = cell2mat(textscan(linecont, '%f')); 
            else
                if params.verbose, fprintf('  NOTICE:  skipping line %d - no output variable declared: %s %s\n', i, lineid, linecont); end;
            end
        case 'NORM'
            if nargout >= 4
                if params.verbose, fprintf('  %s %s\n', lineid, linecont); end;
                info.norm = linecont; 
            else
                if params.verbose, fprintf('  NOTICE:  skipping line %d - no output variable declared: %s %s\n', i, lineid, linecont); end;
            end
        case 'TIDE_SYSTEM'
            if nargout >= 4
                if params.verbose, fprintf('  %s %s\n', lineid, linecont); end;
                info.tide_system = linecont; 
            else
                if params.verbose, fprintf('  NOTICE:  skipping line %d - no output variable declared: %s %s\n', i, lineid, linecont); end;
            end
        case 'END_OF_HEAD'
            header = false;
        otherwise
            if ~header
                if params.verbose, fprintf('  NOTICE:  line %d - unknown id: %s %s\n', i, lineid, linecont); end;
            end
    end
end

%% Third, set empty cells to zero, return coefficients
emptycells = cellfun(@isempty, tmp_coeffs);
tmp_coeffs(emptycells) = {0};

coeffs = cell2mat(tmp_coeffs);  % coefficients in clm-format        

minGOout = minGO;
maxGOout = maxGO;

if params.warning && (sum(emptycells(:)) > 0)
    warning('Some coefficients are not present in the input file. Is this OK? You can suppress this warning by calling the function with parameters "''warning'', false".');
end
if params.verbose, fprintf('  ... finished, minimum read degree/order = %d, maximum read degree/order = %d\n', minGOout, maxGOout); end;


