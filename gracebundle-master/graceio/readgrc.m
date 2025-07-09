function [varargout] = readgrc(mnfname,pth,mn)

% READGRC reads GRACE data and converts the spherical harmonic co-efficient
% data files into readable matrices, i.e., CS-format of the SHBundle
%
% cs            = readgrc(mnfname,pth)
% [cs, meancs]  = readgrc(mnfname,pth,mn)
%
% INPUT
% mnfname   -   Text file with the list of all the names of the GRACE files
%               with separate files for GSM, GAA, GAB, GAC, GAD, &
%               GSM(CALSDV)
% pth      -   Path of the text files if its not the same as the working
%               directory
% mn      -   If you want to calculate the mean type 'mean' else leave it
%               blank
%
% OUTPUT
% cs        -   Cell array of the GRACE Level-2 data [months x 10]
%               [org] [code] [rls] [year] [mnth] [days] [vrsn] [lmax] [cs]
%               [stddev cs]
% meancs    -   Mean of the monthly solutions given in CS-format of the
%               SHBundle
%--------------------------------------------------------------------------

% Author: Balaji Devaraju
% Created on: 8 March 2007, Stuttgart.
%--------------------------------------------------------------------------

tic

if nargin == 1
    if ~ischar(mnfname)
        error('Filename must be a character string. Check input.')
    end
    mn = ' ';
    n = 1;
    pth = [];
elseif nargin == 2
    if ~ischar(mnfname) || ~ischar(pth)
        error('Filename and pathname must be a character string. Check input.')
    end
    mn = ' ';
    n = 1;
elseif nargin == 3
    if ~ischar(mnfname) || ~ischar(pth) || ~ischar(mn)
        error('All inputs must be character strings. Check input.')
    end
    n = 0;
end

if ~isempty(pth) && isempty(regexp(mnfname,pth))
    fid = fopen([pth,mnfname],'r');
else
    fid = fopen(mnfname,'r');
end

fname = {};

itr = true;
while itr
    temp = fgetl(fid);
    if ~ischar(temp) && feof(fid)
        itr = false;
    elseif ischar(temp)
        temp = char(temp);
        if isempty(strfind(temp,pth))
            fexist  = exist([pth,temp],'file');
        else
            fexist = exist(temp,'file');
            [s,e]  = regexp(temp,pth);
            temp   = temp(e+1:end);
        end
        if isequal(fexist,2)
            fname   = [fname; {temp}];
        else
            fprintf('Warning: %s does not exist\n', temp)
        end
    end
end
fclose(fid);

len     = size(fname,1);
grcdata = cell(len,10);
meancs  = 0;

lcnt = unique(round(linspace(1,len,10)));

for i = 1:len
    if ~isempty(find('.'==fname{i}))
        code = [fname{i}(1:3),'SD'];
    else
        code = fname{i}(1:3);
    end
    rls = fname{i}(39:42);
    if str2double(rls) == 1
        days = [str2double(fname{i}(16:18)), str2double(fname{i}(24:26)), str2double(fname{i}(7:10))];
        cnt = round(str2double(fname{i}(24:26))/30);
        if cnt > 0
            mnth = cnt;
            year = str2double(fname{i}(20:23));
        elseif cnt == 0 & str2double(fname{i}(16:18)) >330
            mnth = 12;
            year = str2double(fname{i}(12:15));
        else
            mnth = 1;
            year = str2double(fname{i}(12:15));
        end
    else
        days = [str2double(fname{i}(11:13)), str2double(fname{i}(19:21)), str2double(fname{i}(23:26))];
        if days(3) > 7
            cnt = round(str2double(fname{i}(19:21))/30);
            if cnt > 0
                mnth = cnt;
                year = str2double(fname{i}(7:10));
            elseif cnt == 0 & str2double(fname{i}(11:13))>330
                mnth = 12;
                year = str2double(fname{i}(7:10));
            else
                mnth = 1;
                year = str2double(fname{i}(7:10));
            end
        else
            mnth = round(days(2)/7);
            year = str2double(fname{i}(7:10));
            if mnth < 1
                year = str2double(fname{i}(7:10));
                mnth = 52;
            elseif mnth == 1
                year = str2double(fname{i}(15:18));
            end
        end
    end

    org = fname{i}(28:32);
    if strcmp(org,'EIGEN')
        switch fname{i}(34:37)
        case {'G---'}
            vrsn = 0;
        case {'GK2-'}
            vrsn = 2;
        case {'GOL-'}
            vrsn = 3;
        otherwise
            error('Version unknown')
        end
        %if ~isempty(strfind(fname{i}(34:37),'2'))
        %    vrsn = 2;
        %elseif ~isempty(strfind(fname{i}(34:37),'GM60'))
        %    vrsn = 60;
        %else
        %    vrsn = 0;
        %end
    else
        vrsn = str2double(fname{i}(34:37));
    end
    degree = []; degdot = [];
    order = []; orddot = [];
    clm = []; clmdot = [];
    slm = []; slmdot = [];
    dclm = []; dclmdot = [];
    dslm = []; dslmdot = [];

    fid = fopen([pth,fname{i}]);

    itr = true;
    if isempty(strfind(code,'SD'))
        while itr
            temp = fgets(fid);
            if ~ischar(temp) && feof(fid)
                itr = false;
                break;
            end
            if ~isempty(strfind(temp,'GRCOF2'))
                degree = [degree; str2double(temp(9:11))];
                order = [order; str2double(temp(14:16))];
                clm = [clm; str2double(temp(18:35))];
                slm = [slm; str2double(temp(37:54))];
                dclm = [dclm; str2double(temp(56:65))];
                dslm = [dslm; str2double(temp(67:76))];
            elseif ~isempty(strfind(temp,'GRDOTA'))
                degdot = [degdot; str2double(temp(9:11))];
                orddot = [orddot; str2double(temp(14:16))];
                clmdot = [clmdot; str2double(temp(18:35))];
                slmdot = [slmdot; str2double(temp(37:54))];
                dclmdot = [dclmdot; str2double(temp(56:65))];
                dslmdot = [dslmdot; str2double(temp(67:76))];
            end
        end
        if (min(degree)>0)
            m = find(min(degree)==degree);
            minord = min(order(m));
            if minord ~= 0
                order = [(0:(minord-1))';order];
                degree = [ones(minord-1,1).*min(degree); degree];
                clm = [zeros(minord-1,1); clm];
                slm = [zeros(minord-1,1); slm];
                dclm = [zeros(minord-1,1); dclm];
                dslm = [zeros(minord-1,1); dslm];
            end
            for k = 0:m
                order = [(0:k)';order];
                degree = [ones(k+1,1).*k; degree];
                clm = [zeros(k+1,1); clm];
                slm = [zeros(k+1,1); slm];
                dclm = [zeros(k+1,1); dclm];
                dslm = [zeros(k+1,1); dslm];
            end
        end
        lmax = max(degree);
        cs = sc2cs(clm2sc([degree order clm slm]));
        dcs = sc2cs(clm2sc([degree order dclm dslm]));
        clmcol = length(cs)*orddot+degdot+1;
        slmcol = length(cs)*degdot+orddot+1;
        cs(clmcol) = cs(clmcol) + ((year+mnth/12) - 2000)*clmdot;
        cs(slmcol) = cs(slmcol) + ((year+mnth/12) - 2000)*slmdot;
        dcs(clmcol) = dcs(clmcol) + ((year+mnth/12) - 2000)*dclmdot;
        dcs(slmcol) = dcs(slmcol) + ((year+mnth/12) - 2000)*dslmdot;
        lmax = length(cs) - 1;
    else
        while itr
            temp = fgets(fid);
            if ~ischar(temp) && feof(fid)
                itr = false;
                break;
            end
            if length(temp) > 6
                if strcmp(temp(1:6),'CALSDV')
                    degree = [degree; str2double(temp(9:11))];
                    order = [order; str2double(temp(14:16))];
                    dclm = [dclm; str2double(temp(18:28))];
                    dslm = [dslm; str2double(temp(30:40))];
                end
            end
        end
        if (min(degree)>0)
            m = find(min(degree)==degree);
            minord = min(order(m));
            if minord ~= 0
                order = [(0:(minord-1))';order];
                degree = [ones(minord-1,1).*min(degree); degree];
                dclm = [zeros(minord-1,1); dclm];
                dslm = [zeros(minord-1,1); dslm];
            end
            for k = 0:m
                order = [(0:k)';order];
                degree = [ones(k+1,1).*k; degree];
                dclm = [zeros(k+1,1); dclm];
                dslm = [zeros(k+1,1); dslm];
            end
        end
        lmax = max(degree);
        idx = find((degree==2 & order==1));
        if length(idx) == 2
            degree(idx(2)) = [];
            order(idx(2)) = [];
            dclm(idx(2)) = [];
            dslm(idx(2)) = [];
        end
        dcs = sc2cs(clm2sc([degree order dclm dslm]));
        cs = [];
        lmax = length(dcs) - 1;
    end
    fclose(fid);
    grcdata(i,1:10) = {[org] [code] [rls] [year] [mnth] [days] [vrsn] [lmax] [cs] [dcs]};

    if strcmp(mn,'mean')
        if days(3)<40
            if ~strcmp(org,'EIGEN')
                if meancs == 0
                    meancs = meancs+cs;
                elseif length(meancs) < length(cs)
                    meancs = meancs + cs(1:length(meancs),1:length(meancs));
                else
                    meancs = meancs(1:(lmax+1),1:(lmax+1)) + cs;
                end
                n = n+1;
            elseif strcmp(org,'EIGEN')
                if ((year+mnth/12)<(2004+7/12)) | ((year+mnth/12)>(2004+10/12))
                    if meancs == 0
                        meancs = meancs+cs;
                    elseif length(meancs) < length(cs)
                        meancs = meancs + cs(1:length(meancs),1:length(meancs));
                    else
                        meancs = meancs(1:(lmax+1),1:(lmax+1)) + cs;
                    end
                    n = n+1;
                end
            end
        else
            fprintf('Found a geoid \n')
        end
    end

    if ismember(i,lcnt)
        fprintf('.')
    end
end

varargout{1} = grcdata;
if nargout == 2
    varargout{2} = meancs./n;
end

fprintf(' %s \t %g[s]\n', code, toc)
