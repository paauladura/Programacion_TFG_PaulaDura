function [cs, ncs, header, cst, ncst] = readicgem(filename)
    
% READICGEM reads ICGEM format files containing spherical harmonic coefficients
% 
% [shcoeff, header] = readicgem(filename)
%
% INPUT
% filename - Filenmae as a string
%
% OUTPUT
% cs, ncs   - Spherical harmonic coefficients of the mean field
% header    - A structure file containing the header information available in
%             the given file.
% cst, ncst - Spherical harmonic coefiicients of the trend estimates
%-------------------------------------------------------------------------------
%

% Author: Balaji Devaraju
% Date: 9 October 2021
% Place: IIT Kanpur, India
%-------------------------------------------------------------------------------

narginchk(1, 1)
nargoutchk(1, 5)

if ischar(filename)
    if ~isequal(exist(filename), 2)
        error(sprintf('%s does not exist', filename))
    end
else
    error('Filename is not a string')
end

fid = fopen(filename, 'r+');
tic
C = cell(1e5,3);
n = 0;
begin_of_head   = false;
end_of_head     = false;
header = [];
while 1
    one_line = fgets(fid);
    
    if ~ischar(one_line) && feof(fid)
        break
    end
    
    cells = strsplit(one_line);
    key   = lower(cells{1});
    
    if strcmp(key, '')
        continue
    end
    
    if end_of_head
        n = n + 1;
        
        cells = cells(~strcmp('', cells));
        
        C{n,1} = cells{1};
        if hasErrors
            C{n,2} = cellfun(@str2double, cells(2:7));
            C{n,3} = cellfun(@str2double, cells(8:end));
        else
            C{n,2} = cellfun(@str2double, cells(2:5));
            C{n,3} = cellfun(@str2double, cells(6:end));
        end
        
    elseif begin_of_head && ~end_of_head
        switch key
            case {'max_degree', 'maxdegree'}
                header.max_degree = str2double(cells{2});
            case {'radius'}
                header.ae = str2double(cells{2});
            case {'earth_gravity_constant'}
                header.GM = str2double(cells{2});
            case {'modelname', 'model_name'}
                header.modelname = strjoin(cells(2:end));
            case { 'generating_institute' }
                header.institute = strjoin(cells(2:end));
            case {'norm', 'normalization', 'normalisation'}
                header.norm = strjoin(cells(2:end));
            case {'errors'}
                if(strcmp(cells{2}, 'no'))
                    hasErrors = false;
                else
                    hasErrors = true;
                end
                header.errors = strjoin(cells(2:end));
            case {'tide_system', 'tidesystem'}
                header.tide_system = strjoin(cells(2:end));
            case {'time_period_of_data', 'time_period', 'timeperiod'}
                header.time_period = strjoin(cells(2:end));
            case {'ref_epoch_[mjd]'}
                header.ref_epoch_mjd = strjoin(cells(2:end));
            case {'ref_epoch_[yyyy/mm/dd/h]'}
                header.ref_epoch_cal = strjoin(cells(2:end));
            case {'key'}
                cells       = cells(~strcmp('', cells));
                header.key  = cells(2:end);
            case {'end_of_head'}
                end_of_head = true;
            otherwise
                header.(key) = strjoin(cells(2:end));
        end
        
    else
        idx     = ~strcmp('', cells);
        cells   = cells(idx);
        
        switch key
            case {'begin_of_head'}
                begin_of_head = true;
            case {'time_coverage_start', 'time_coverage_end'}                
                header.(key) = cells{3};
                if length(cells) >= 4
                    header.(key) = [header.(key), 'T', cells{4}];
                end
            case { 'generating_institute' }
                header.institute = strjoin(cells(2:end));
                
            case {'time_period_of_data', 'time_period', 'timeperiod'}
                header.time_period = horzcat(cells{2:end});
            case {'ref_epoch_[mjd]'}
                header.ref_epoch_mjd = cells{2};
            case {'ref_epoch_[yyyy/mm/dd/h]'}
                header.ref_epoch_cal = cells{2};
        end
    end
end
fprintf('%s read in %g seconds\n', filename, toc);
fclose(fid);

C       = C(1:n, :);
codes   = unique(C(1:n, 1));

for k = 1:length(codes)
    idx = strcmp(codes{k}, C(:,1));
    tmp = cell2mat(C(idx,2));
    cs  = sc2cs(clm2sc(tmp(:,1:4)));
    dcs = sc2cs(clm2sc(tmp(:,[1:2, 5:6])));
    
    shcoeff.(codes{k}).cs       = cs;
    if hasErrors
        shcoeff.(codes{k}).errors   = dcs;
    end
    m = find(idx);
    if ~isnan(C{m(k),3}) & ~strcmp(C{m(k),3}, '')
        shcoeff.(codes{k}).info = [tmp(:,1:2), cell2mat(C(idx,3))];
    end
end

if any(strcmp('gfct', codes))
    cs  = blkdiag(shcoeff.gfc.cs, zeros(size(shcoeff.gfct.cs,1)-size(shcoeff.gfc.cs,1))) + shcoeff.gfct.cs;
    ncs = shcoeff.gfct.errors;
else
    cs  = shcoeff.gfc.cs;
    ncs = shcoeff.gfc.errors;
end

if any(strcmp('trnd', codes))
    cst = shcoeff.trnd.cs;
    ncst= shcoeff.trnd.errors;
else
    cst = [];
    ncst= [];
end