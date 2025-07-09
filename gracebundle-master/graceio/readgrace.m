function [gxx_cell, GM, ae] = readgrace(pth, gxx)

% READGRACE reads GRACE spherical harmonic data as provided by the three
% data processing centers stored in a directory/folder
%
% gxx_cell = readgrace(pth, gxx)
%
% INPUT
% pth 	- Path where all the GRACE spherical harmonic data is stored
% gxx 	- Character code for reading specific files
% 			'GAA' - Non-tidal atmosphere
% 	 		'GAB' - Non-tidal ocean
% 			'GAC' - 'GAA'+'GAB'
% 			'GAD' - Ocean bottom pressure
% 			'GSM' - GRACE monthly gravity field solutions
% 			'GSD' - GRACE calibrated errors
%
% OUTPUT
% gxx_cell - Cell array of the GRACE Level-2 data [months x 10]
%               [org] [code] [rls] [year] [mnth] [days] [vrsn] [lmax] 
% 				[cs] [stddev cs]
%-----------------------------------------------------------------------

% Check number of input arguments
narginchk(2, 2);

% Check inputs
if ischar(pth) && ischar(gxx)
	if isequal(exist(pth), 7)
		if ~strcmp(pth(end), '/')
			pth = [pth, '/'];
		end
	else
		err = sprintf('%s directory does not exist', pth);
		error(err)
	end
else
	error('Inputs must be character strings. Verify input!')
end

% List all the files
switch gxx
	case { 'GAA', 'GAB', 'GAC', 'GAD' }
		list_ = dir([pth, gxx, '*05.gz']);
	case { 'GSM' }
		list_ = dir([pth, gxx, '*5.gz']);
		if isempty(list_)
			list_ = dir([pth, gxx, '*5a.gz']);
		end
	case { 'GSD' }
		list_ = dir([pth, 'GSM*5.txt']);
		if isempty(list_)
			list_ = dir([pth, 'GSM*5a.txt']);
		end
	otherwise
		err = sprintf('Unknown choice "%s"', gxx);
end

% Check if there are any files in the given directory
if isempty(list_)
	error('No files found! Please verify your inputs.')
else
	list_ = cellstr(vertcat(list_(:).name));
    list_ = strcat(pth, list_);
end

[gxx_cell, GM, ae] = cellfun(@extract_data, list_, 'UniformOutput', false);

gxx_cell = vertcat(gxx_cell{:});
GM       = GM{1};
ae       = ae{1};

end


% Internal function
function [gxx, GM, ae] = extract_data(fnm)

fprintf('%s\n', fnm);

% Current working directory
id = regexp(fnm, '/');
curr_dir = fnm(1:id(end));

if isequal(exist([curr_dir,'gxx']), 7)
    curr_dir = [curr_dir, 'gxx/'];
else
    curr_dir = [curr_dir, 'gxx/'];
    mkdir(curr_dir);
end

% Unzip the file
if ~isempty(regexp(fnm, '.gz'))
    copyfile(fnm, curr_dir);
    fnm = [curr_dir, fnm(id(end)+1:end)];
    fnm = gunzip(fnm);
    fnm = fnm{1};
end

% Find the underscores in the filename
s = regexp(fnm, '_');
e = s-1;
e = [e, regexp(fnm, '\.')-1];
s = [1, s+1];

% Data type GXX
id = regexp(fnm, '/');
if ~isempty(regexp(fnm, '.txt'))
    code = [fnm(id(end)+1:id(end)+3), 'SD'];
else
    code = fnm(id(end)+1:id(end)+3);
end

% Time span
tmp     = fnm(s(2):e(2));
idx     = regexp(tmp, '-');
yr      = str2double(tmp(1:4));
days    = [str2double(tmp(5:7)), ...
            str2double(tmp(13:15)), ...
            str2double(fnm(s(3):e(3)))];
mnth    = doy2cal(yr, days(1));

% Processing centre
org     = fnm(s(4):e(4));
if strcmp(org,'EIGEN')
    switch fnm(s(5):e(5))
        case {'G---'}
            vrsn = 0;
        case {'GK2-'}
            vrsn = 2;
        case {'GOL-'}
            vrsn = 3;
        otherwise
            error('Version unknown')
    end
else
    vrsn = 0;
end

% Release number
rls = fnm(s(end):e(end));
s = regexp(rls, '0');
rls = rls(s(end)+1:end);

% Read the file
fid = fopen(fnm, 'r+');

loop = true;
clm  = cell(sum(1:2301), 1);
dclm = clm;

GM = NaN;
ae = NaN;
lmax = NaN;

k = 0;
while loop
    tmp = fgetl(fid);
    k   = k + 1;

    if ischar(tmp) && ~feof(fid)
        
        if ~isempty(tmp)
            str = textscan(tmp, '%s');
        else
            str{1} = {};
        end

        if ~isempty(str{1})
            switch str{1}{1}
            case { 'EARTH' }
                GM = str2double(str{1}{2});
                ae = str2double(str{1}{3});
            case { 'SHM' }
                lmax = str2double(str{1}{2});
            case { 'GRCOF2' }
                clm{k} = cellfun(@str2double, str{1}(2:5)');
                dclm{k} = cellfun(@str2double, str{1}(6:7)');
            case { 'CALSDV' }
                dclm{k} = cellfun(@str2double, str{1}(2:5)');
                clm{k} = dclm{k};
                clm{k}(3:4) = 0;
            case { 'FIRST', 'CMMNT' }
                fprintf('%s\n', tmp);
            end
        end
    else
        loop = false;
    end
end

% Reformat the cell into a matrix
idx = cellfun(@isempty, clm);
clm(idx) = [];
clm = vertcat(clm{:});

idx = cellfun(@isempty, dclm);
dclm(idx) = [];
dclm = vertcat(dclm{:});

% Whne reading calibrated standard deviations it is not possible to find
% out the maximum degree of spherical harmonic expansion while reading
% the file. We do it after reading all the data from the file.
if ~isempty(regexp(fnm, '.txt'))
    lmax = max(dclm(:,1));
end

% Convert [l m Clm Slm] to [C\S]
dclm = sc2cs(clm2sc([clm(:,1:2) dclm]));
clm  = sc2cs(clm2sc(clm));

gxx = {org, code, rls, yr, mnth, days, vrsn, lmax, clm, dclm};
end
