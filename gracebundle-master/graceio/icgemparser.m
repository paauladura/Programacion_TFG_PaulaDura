function [grav] = icgemparser(filename, varargin)

% ICGEMPARSER reads files in the ICGEM format.
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
%   grav ............ Gravity field model as a structure file
% 						max_degree  Maximum degree of spherical 
%                                   harmonic expansion
% 						ae          Semi-major axis of the reference 
%                                   earth ellipsoid
% 						GM          Gravitational constant
% 						modelname   Name of the model
% 						norm        Normalization used
% 						tide_system Tide system 'zero tides', 'tide free'
% 						gfc         Gravity field coefficients in C\S 
%                                   format
% 					        knm     Estimates
% 						    dknm    Errors
%
% EXAMPLE:
%   icgemparser('./GRACE/ITG-Grace2010s.gfc'); 
%   icgemparser('./GRACE/ITG-Grace2010s.gfc', 'max_lm', 20); 
%   icgemparser('./GRACE/ITG-Grace2010s.gfc', 'max_lm', 20, 'min_lm', 5, 'verbose', true); 
%
% SCREEN DUMP:
%    'comment' lines, as well as unknown line identifiers
%
% USES:
%    uberall/GETOPT

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

% read header
fid       = fopen(filename, 'r');
hasErrors = true;
hlines    = 0;
while 1
    one_line = fgets(fid);
    if ~ischar(one_line), break, end
    hlines  = hlines + 1;
    if isempty(one_line)
        continue
    else
        keyword = textscan(one_line,'%s',1);
        if isempty(keyword{1})
            continue
        else
            keyword = lower(keyword{1}{1});
        end
    end


    switch keyword
    case {'max_degree', 'maxdegree'}
        cells           = textscan(one_line,'%s%f',1);
        grav.max_degree = cells{2};
    case {'radius'}
        cells   = textscan(one_line,'%s%f',1);
        grav.ae = cells{2};
    case {'earth_gravity_constant'}
        cells   = textscan(one_line,'%s%f',1);
        grav.GM = cells{2};
    case {'modelname', 'model_name'}
        cells          = textscan(one_line,'%s%s',1);
        grav.modelname = cells{2}{1};
    case { 'generating_institute' }
        cells = strread(one_line, '%s');
        grav.institute = horzcat(cells{2:end});
    case {'norm', 'normalization', 'normalisation'}
        cells     = textscan(one_line,'%s%s',1);
        grav.norm = cells{2}{1};
    case {'errors'}
        cells = textscan(one_line,'%s%s',1);
        if(strcmp(cells{2}{1}, 'no'))
          hasErrors = false;
        end
    case {'tide_system', 'tidesystem'}
        cells = textscan(one_line,'%s%s',1);
        grav.tide_system = cells{2}{1};
    case {'time_period_of_data', 'time_period', 'timeperiod'}
        cells = textscan(one_line, '%s%s%s%s', 1);
        grav.time_period = vertcat(cells{2:end});
    case {'ref_epoch_[mjd]'}
        cells = textscan(one_line, '%s%f', 1);
        grav.ref_epoch_mjd = cells{2};
    case {'ref_epoch_[yyyy/mm/dd/h]'}
        cells = textscan(one_line, '%s%s', 1);
        grav.ref_epoch_cal = cells{2}{1};
    case {'end_of_head', 'end of head'}
        fprintf('Found %s at line #%g\n', keyword, hlines);
        break
    end
end
fclose(fid);

lmin = params.min_lm + 1;
if isinf(params.max_lm)
    lmax = grav.max_degree + 1;
else
    lmax = params.max_lm + 1;
end

if hasErrors
    
    if isoctave
        [code, l, m, C, S, dC, dS] = textread(filename, '%s %f %f %f %f %f %f', 'headerlines', hlines);
    else
        [code, l, m, C, S, dC, dS] = textread(filename, '%s %f %f %f %f %f %f %*[^\n]', 'headerlines', hlines);
    end
    
    codes = unique(code);
    
    for k = 1:length(codes)
        idx  = ismember(code, codes{k});
        mtmp = m(idx);
        ltmp = l(idx);
        midx = (mtmp > 0);
        tmp  = S(idx);
        grav.(codes{k}).knm = zeros(grav.max_degree+1);
        grav.(codes{k}).knm((l(idx) + 1) + (m(idx) * (grav.max_degree + 1))) = C(idx);
        grav.(codes{k}).knm(mtmp(midx) + (ltmp(midx)*(grav.max_degree + 1))) = tmp(midx);
        grav.(codes{k}).knm = grav.(codes{k}).knm(lmin:lmax,lmin:lmax);

        tmp  = dS(idx);
        grav.(codes{k}).dknm = zeros(grav.max_degree+1);
        grav.(codes{k}).dknm((l(idx) + 1) + (m(idx) * (grav.max_degree + 1))) = dC(idx);
        grav.(codes{k}).dknm(mtmp(midx) + (ltmp(midx)*(grav.max_degree + 1))) = tmp(midx);
        grav.(codes{k}).dknm = grav.(codes{k}).dknm(lmin:lmax,lmin:lmax);
    end
else
    if isoctave
        [code, l, m, C, S] = textread(filename, '%s %f %f %f %f', 'headerlines', hlines);
    else
        [code, l, m, C, S] = textread(filename, '%s %f %f %f %f %*[^\n]', 'headerlines', hlines);
    end
    
    codes = unique(code);
    for k = 1:length(codes)
        idx  = ismember(code, codes{k});
        mtmp = m(idx);
        ltmp = l(idx);
        midx = (mtmp > 0);
        tmp  = S(idx);
        grav.(codes{k}).knm = zeros(grav.max_degree+1);
        grav.(codes{k}).knm((l(idx) + 1) + (m(idx) * (grav.max_degree + 1))) = C(idx);
        grav.(codes{k}).knm(mtmp(midx) + (ltmp(midx)*(grav.max_degree + 1))) = tmp(midx);
        grav.(codes{k}).knm = grav.(codes{k}).knm(lmin:lmax,lmin:lmax);
    end
end
