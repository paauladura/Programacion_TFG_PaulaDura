function [fvar, thetaRAD, lamRAD, varargout] = gshscov(field,lmax,quant,grid,n,fil,h,ivartype,ovartype)

% GSHSCOV propagates the covariance information of a spherical harmonic
% field to the spatial domain based on the GSHS approach.
% The methods are based on the 2-D Fourier approach as shown by R. Haagmans
% and M. van Gelderen, 1991. JGR (B).
%
% [fvar,thetaRAD,lamRAD] = gshscov(field,lmax,quant,grid,n,fil,h,ivartype,ovartype)
%
% IN:
%    field .... A [(L+1)^2 * (L+1)^2] variance-covariance matrix arranged
%               in Colombo ordering, or, a structure with two variables:
%               1. a character array with the covariance matrix file;
%               2. path where the file is stored;
%               or, a character array of the filename and the path
%               together,
%               or, if only the variances are available then provide them
%               in either SC, CS, or Colombo ordering formats, and provide
%               an empty matrix for ivartype.
%               Use the first option only when the maximum degree of the
%               spherical harmonic expansion is [lmax <= 70].
%    lmax ..... Maximum degree of the spherical harmonic expansion.
%    quant .... optional string argument, defining the field quantity:
%               'none' (coefficients define the output), 'geoid',
%               'potential', 'dg' or 'gravity' (grav. anomaly), 'tr' (grav.
%               disturbance), or 'trr' (2nd rad. derivative), 'water'
%               (water equivalent heights), 'smd' (surface mass density).
%               If in the case of cross-covariance computation quant is a
%               cell array with two-variables. One for each field.
%                                                           -def: 'none'
%    grid ..... optional string argument, defining the grid:
%               1. 'pole' or 'mesh': equi-angular (N+1)*2N, including poles
%               and Greenwich meridian.                     -def: 'pole'
%               2. 'block' or 'cell': equi-angular block mid-points.
%               3. 'neumann' or 'gauss': Gauss-grid (N+1)*2N
%    n ........ grid size parameter N. #longitude samples: 2N
%               #latitude samples N ('blocks') or N+1.		- def: N=L
%    fil ...... filter types that are used for spatial averaging. Input is
%               a structure variable with the name and the required
%               parameters of the filters.
%               options:
%               1. 'gauss'  - Gaussian smoothing. Parameter for the filter
%               is a smoothing radius. [radians]
%               2. 'pell'   - Pellinen smoothing. Parameter for the filter
%               is a spherical cap angle psi. 0 < psi < pi [radians]
%               3. 'han'    - Non-isotropic smoothing proposed by Shin-Chan
%               Han et al., 2005, GJI 163, 18--25. Parameters for the
%               filter are
%                       m1  - order threshold [no units]
%                       r0  - fundamental radius [km]
%                       r1  - secondary radius [km]
%               4. 'other'   -  Other filter co-efficients can either be 
%               given as a [l Wl] vector (isotropic case) or as a matrix 
%               in CS- / SC- / [l m Clm Slm] in colombo ordering- formats
%               (anisotropic case).
%               5. 'hann'   - Hanning smoothing. Parameter is smoothing
%               radius r [radians]. The radius defined here is different from
%               Gaussian, because at the radius specified the filter
%               co-efficients become zero. Whereas in Gaussian smoothing
%               at the radius the filter coefficient value becomes half.
%               6. 'none'   - no filtering will be applied
%                                                           - def: 'none'
%    h ........ height for upward continuation [m]. It is also possible to
%               provide to different heights if the covariance matrix is
%               between two fields that are representative of two different
%               heights.
%                                                           - def: 0.
%    ivartype.. type of variance-covariance propagation
%               1. 'full' - full spectral variance-covariance matrix is
%                           propagated.
%               2. 'block'- block-diagonal variance-covariance matrix is
%                           propagated.
%               3. 'diag' - only the diagonal of the variance-covariance
%                           matrix is propagated.
%                                                           - def: 'diag'
%    ovartype.. type of propagated variance-covariance information
%               1. 'full' - varainces and covariances of all the points are
%               calculated and stored in *.mat files. The files contain the
%               inter-latitude covariances.
%               2. 'diag' - only variances are stored. The final output is
%               a matrix with only variances for each grid point.
%
% OUT:
%    fvar ..... the global variances of the field
%    thetaRAD . vector of co-latitudes [radians]
%    lamRAD ... vector of longitudes [radians]
%    cvflist .. a text file with the list of all the covariance matrix
%               files.
%
% USES:
%   SHbundle\cs2sc, 
%               cssc2clm,
%               covord,
%               eigengrav, 
%               plm, 
%               clm2sc, 
%               sc2cs
%   FilterBundle\pellinen, 
%               gaussfltr, 
%               vonhann,
%               hannoniso
%   uberall\isint,
%           standing,
%           grule
%
% SEE ALSO:
%    gshs_, gsha

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2021-04-12: MA, remove revision statement on deprecated function 
%    2014-07-24: BD, Filter radii inputs only in radians
%    2007-11-28: BD, initial version
%-----------------------------------------------------------------------------
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
% ----------------------------------------------------------------------------

tic
fprintf('------------------------------------------------------------------\n')
fprintf('\t\t Propagating spectral covariances to space \n')
fprintf('------------------------------------------------------------------\n')

fprintf('Checking input arguments ... ')

%-----------------------------------------------------
% Checking the input and augmenting initial values
%-----------------------------------------------------
if nargin > 9, error('Too many input arguments!'),     end
if nargin < 9 || isempty(ovartype), ovartype = 1;      end
if nargin < 8,                      ivartype = 2;      end
if nargin < 7 || isempty(h),        h        = 0;      end
if nargin < 6 || isempty(fil),      fil      = struct('type','none'); end
if nargin < 5 || isempty(n),        n        = lmax;   end
if nargin < 4 || isempty(grid),     grid     = 'mesh'; end
if nargin < 3 || isempty(quant),    quant    = 'none'; end
if nargin < 2, error('Insufficient input arguments'),  end

if ~(isscalar(n) && isint(n)), error('n must be integer & scalar'), end
if ~ischar(grid), error('grid argument must be string'), end
grid  = lower(grid);

%-------------------------
% Checking field types
%-------------------------
if isstruct(field)
    fn          = fieldnames(field);
    fname       = field.(fn{1});
    pth         = field.(fn{2});
    delfile     = 0;
elseif ischar(field)
    fname       = field;
    pth         = pwd;
    pth         = [pth,'/'];
    delfile     = 0;
elseif isempty(ivartype)
    ivartype = 'diag';
    [rows,cols] = size(field);
    if isequal(cols,4) && isequal(rows, sum(1:lmax+1))
        field = sc2cs(clm2sc(field));
    elseif isequal(cols,(2*lmax+1))
        field = sc2cs(field);
    elseif ~isequal(rows,cols)
        error('Error matrix not in the SC, CS, or Colombo ordering formats')
    end
    vec = cssc2clm(field,lmax);
    vec = sparse(diag([vec(:,3); vec(lmax+2:end,4)]));
    save('qmatnow.mat', 'vec');
    % vec         = []; 
    fname       = 'qmatnow.mat';
    pth         = pwd;
    pth         = [pth, '/'];
    delfile     = 1;
else
    [rows, cols] = size(field);
    if ~isequal(rows, cols) || ~isequal((lmax+1)^2, rows)
        error('Matrix dimensions do not agree with the maximum degree of the spherical harmonic expansion')
    end
    save 'qmatnow.mat' field
    fname       = 'qmatnow.mat';
    pth         = pwd;
    pth         = [pth,'/'];
    delfile     = 1;
end

%----------------------------------------------------
% Checking required input and output matrix types
%----------------------------------------------------
if nargin == 8 || nargin == 9
    if strcmp(ivartype,'full')
        ivartype = 0;
    elseif strcmp(ivartype,'block')
        ivartype = 1;
    elseif strcmp(ivartype,'diag')
        ivartype = 2;
    else
        display('Warning: Unknown covariance matrix structure')
    end
end

if nargin == 9
    if strcmp(ovartype,'full')
        ovartype = 0;
        mkdir(pth,'covmats')
        covpth = [pth, 'covmats/'];
    elseif strcmp(ovartype,'diag')
        ovartype = 1;
    else
        display('Warning: Unknown covariance matrix structure')
    end
end

fprintf('done\t')
toc

tic
fprintf('Preparing for spherical harmonic synthesis ... ')
% ------------------------
% Grid definition.
% ------------------------
dt = 180/n;
if strcmp(grid,'pole') || strcmp(grid,'mesh')
    theta = (0:dt:180)';
    lam   = (0:dt:360-dt);
elseif strcmp(grid,'block') || strcmp(grid,'cell')
    theta = (dt/2:dt:180)';
    lam   = (dt/2:dt:360);
elseif strcmp(grid,'neumann') || strcmp(grid,'gauss')
    gx = grule(n+1);
    theta = flipud(acos(standing(gx)))*180/pi;
    lam   = (0:dt:360-dt);
else
    error('Unknown grid type')
end
nlat = length(theta);
nlon = length(lam);

%------------------------------------
% Preprocessing on the coefficients:
%    - specific transfer (quant)
%    - filtering (fil,cap)
%    - upward continuation (h)
% -----------------------------------

l    = (0:lmax)';
transf = eigengrav(l,quant,h);

%-----------------------------------
% Processing filter (coefficients)
%-----------------------------------

fldvar = fieldnames(fil);
if strcmp(fil.(fldvar{1}),'gauss')
    filcff = gaussfltr(fil.(fldvar{2}),lmax);
    filcff = filcff.*transf;
elseif strcmp(fil.(fldvar{1}),'pell')
    filcff = pellinen(l,fil.(fldvar{2}));
    filcff = filcff.*transf;
elseif strcmp(fil.(fldvar{1}),'hann')
    filcff = vonhann(fil.(fldvar{2}),lmax);
    filcff = filcff.*transf;
elseif strcmp(fil.(fldvar{1}),'han')
    m1      = fil.(fldvar{2});
    r0      = fil.(fldvar{3});
    r1      = fil.(fldvar{4});
    filcff = sc2cs(clm2sc(hannoniso(r0,m1,r1,lmax)));
    filcff = filcff.*sc2cs(clm2sc(cssc2clm([(0:lmax)' transf],lmax)));
elseif strcmp(fil.(fldvar{1}),'other')
    filcff = fil.(fldvar{2});
    filcff = sc2cs(clm2sc(cssc2clm(filcff,lmax)));
    filcff = filcff.*sc2cs(clm2sc(cssc2clm([(0:lmax)' transf],lmax)));
elseif strcmp(fil.(fldvar{1}),'none')
    filcff = transf;
end

%------------------------------------------
% Calculating the cc cs sc ss coefficients
%------------------------------------------
cols = lmax + 1;

% If N > L this will appear as a 2D zero-padding as we are dealing with 2D
% Fourier transforms
cc = zeros(cols);
cs = zeros(cols);
sc = zeros(cols);
ss = zeros(cols);

fprintf('done\t')
toc

tic
fprintf('Starting spherical harmonic synthesis of ')

% Initializing variance matrix
fvar = zeros(nlat,nlon);

hwb    = waitbar(0,'Percentage of latitude points ready ...');
set(hwb,'NumberTitle','off','Name','GSHS - Covariance propagation')

trad = theta * pi/180;
% Writing covariance files
if ivartype==0
	fprintf('full covariance matrix ... ')
    covcell = covord(fname,pth,lmax,0);
    if delfile == 1, delete(fname), end

    if size(filcff,2) == 1
        for k = 1:lmax+1
            for m = k:lmax+1
                temp = (filcff(m:end)*filcff(k:end)');
                covcell.cc{m,k} = temp.*covcell.cc{m,k};
                covcell.ss{m,k} = temp.*covcell.ss{m,k};
                covcell.sc{m,k} = temp.*covcell.sc{m,k};
                covcell.cs{m,k} = temp.*covcell.cs{m,k};
                if m~=k
                    temp = (filcff(k:end)*filcff(m:end)');
                    covcell.cc{k,m} = temp.*covcell.cc{k,m};
                    covcell.ss{k,m} = temp.*covcell.ss{k,m};
                    covcell.sc{k,m} = temp.*covcell.sc{k,m};
                    covcell.cs{k,m} = temp.*covcell.cs{k,m};
                end
            end
        end
    elseif size(filcff,1)==size(filcff,2) && size(filcff,1)==lmax+1
        for k = 1:lmax+1
            for m = k:lmax+1
                temp = (filcff(m:end,m)*filcff(k:end,k)');
                covcell.cc{m,k} = temp.*covcell.cc{m,k};
                covcell.ss{m,k} = temp.*covcell.ss{m,k};
                covcell.sc{m,k} = temp.*covcell.sc{m,k};
                covcell.cs{m,k} = temp.*covcell.cs{m,k};
                if m~=k
                    temp = (filcff(k:end,k)*filcff(m:end,m)');
                    covcell.cc{k,m} = temp.*covcell.cc{k,m};
                    covcell.ss{k,m} = temp.*covcell.ss{k,m};
                    covcell.sc{k,m} = temp.*covcell.sc{k,m};
                    covcell.cs{k,m} = temp.*covcell.cs{k,m};
                end
            end
        end
    end
    if ovartype==0,
        fid = fopen([covpth,'cvflist.txt'],'W+');
        % Calculating the 2-D Fourier coefficients
        %
        for i = 1:nlat
            for j = i:nlat
                for k = 0:lmax
                    for m = k:lmax
                        if i ~= j
                            pnmi = plm(m:lmax,m,trad([i;j]));
                            pnmj = pnmi(2,:);
                            pnmi(2,:) = [];
                            if m~=k
                                plki = plm(k:lmax,k,trad([i;j]));
                                plkj = plki(2,:);
                                plki(2,:) = [];
                            elseif m == k
                                plki = pnmi;
                                plkj = pnmj;
                            end
                        elseif i == j
                            pnmi = plm(m:lmax,m,trad(i));
                            pnmj = pnmi;
                            if m~=k
                                plki = plm(k:lmax,k,trad(i));
                                plkj = plki;
                            elseif m == k
                                plki = pnmi;
                                plkj = pnmj;
                            end
                        end
                        temp = pnmi'.* plkj;
                        cc(m+1,k+1) = sum(sum(temp.* covcell.cc{m+1,k+1}));
                        cs(m+1,k+1) = sum(sum(temp.* covcell.cs{m+1,k+1}));
                        sc(m+1,k+1) = sum(sum(temp.* covcell.sc{m+1,k+1}));
                        ss(m+1,k+1) = sum(sum(temp.* covcell.ss{m+1,k+1}));

                        if m~=k
                            temp = plki'*pnmj;
                            cc(k+1,m+1) = sum(sum(temp.* covcell.cc{k+1,m+1}));
                            cs(k+1,m+1) = sum(sum(temp.* covcell.cs{k+1,m+1}));
                            sc(k+1,m+1) = sum(sum(temp.* covcell.sc{k+1,m+1}));
                            ss(k+1,m+1) = sum(sum(temp.* covcell.ss{k+1,m+1}));
                        end
                    end
                end
                mlam = (lam'*pi/180)*(0:cols-1);
                f = cos(mlam)*cc*cos(mlam');
                f = f + sin(mlam)*sc*cos(mlam');
                f = f + cos(mlam)*cs*sin(mlam');
                f = f + sin(mlam)*ss*sin(mlam');
                % Saving the covariance fields in files
                if floor(theta(i)) < 10,
                    lat1 = ['00',num2str(theta(i),'%7.3f')];
                elseif floor(theta(i)) >= 10 && floor(theta(i)) < 100
                    lat1 = ['0',num2str(theta(i),'%7.3f')];
                elseif floor(theta(i)) >= 100
                    lat1 = num2str(theta(i),'%7.3f');
                end
                if i ~= j
                    if floor(theta(j)) < 10,
                        lat2 = ['00',num2str(theta(j),'%7.3f')];
                    elseif floor(theta(j)) >= 10 && floor(theta(j)) < 100
                        lat2 = ['0',num2str(theta(j),'%7.3f')];
                    elseif floor(theta(j)) >= 100
                        lat2 = num2str(theta(j),'%7.3f');
                    end
                elseif i == j
                    lat2 = lat1;
                    fvar(i,:) = diag(f)';
                end
                save([covpath, 'cov', lat1,'-', lat2, '.mat'], 'f')
                fprintf(fid,'%s \n',['cov', lat1,'-', lat2, '.mat']);
                cc = zeros(cols);
                cs = zeros(cols);
                sc = zeros(cols);
                ss = zeros(cols);
            end
            waitbar(i/nlat)
        end
        fclose(fid);
        fprintf('done\t')
    elseif ovartype == 1
        for i = 1:nlat
            for k = 0:lmax
                for m = k:lmax
                    pnm = plm(m:lmax,m,trad(i));
                    if m~=k
                        plk = plm(k:lmax,k,trad(i));
                    elseif m == k
                        plk = pnm;
                    end
                    temp = pnm'*plk;
                    cc(m+1,k+1) = sum(sum(temp.* covcell.cc{m+1,k+1}));
                    cs(m+1,k+1) = sum(sum(temp.* covcell.sc{k+1,m+1}'));
                    sc(m+1,k+1) = sum(sum(temp.* covcell.sc{m+1,k+1}));
                    ss(m+1,k+1) = sum(sum(temp.* covcell.ss{m+1,k+1}));

                    if m~=k
                        cc(k+1,m+1) = cc(m+1,k+1);
                        cs(k+1,m+1) = sc(m+1,k+1);
                        sc(k+1,m+1) = cs(m+1,k+1);
                        ss(k+1,m+1) = ss(m+1,k+1);
                    end
                end
            end
            mlam = (lam'*pi/180)*(0:cols-1);
            f = cos(mlam)*cc*cos(mlam');
            f = f + sin(mlam)*sc*cos(mlam');
            f = f + cos(mlam)*cs*sin(mlam');
            f = f + sin(mlam)*ss*sin(mlam');
            fvar(i,:) = diag(f)';
            cc = zeros(cols);
            cs = zeros(cols);
            sc = zeros(cols);
            ss = zeros(cols);
            waitbar(i/nlat)
        end
        fprintf('done\t')
    end
elseif ivartype==1
	fprintf('block-diagonal covariance matrix ... ')
    covcell = covord(fname,pth,lmax,1);
    if delfile == 1, delete(fname), end

    if size(filcff,2) == 1
        for m = 1:lmax+1
            temp = filcff(m:end);
            temp = temp*temp';
            covcell.cc{m} = temp.*covcell.cc{m};
            covcell.ss{m} = temp.*covcell.ss{m};
        end
    elseif size(filcff,1)==size(filcff,2) && size(filcff,1)==lmax+1
        for m = 1:lmax+1
            temp = filcff(m:end,m);
            temp = temp*temp';
            covcell.cc{m} = temp.*covcell.cc{m};
            covcell.ss{m} = temp.*covcell.ss{m};
        end
    end
    if ovartype == 0
        fid = fopen([covpth,'cvflist.txt'],'W+');
        % Calculating the 2-D Fourier coefficients
        %
        for i = 1:nlat
            for j = i:nlat
                for m = 0:lmax
                    if i ~= j
                        pnmi = plm(m:lmax,m,trad([i;j]));
                        pnmj = pnmi(2,:);
                        pnmi(2,:) = [];
                    elseif i == j
                        pnmi = plm(m:lmax,m,trad(i));
                        pnmj = pnmi;
                    end
                    temp = pnmi'*pnmj;
                    cc(m+1,m+1) = sum(sum(temp.* covcell.cc{m+1}));
                    ss(m+1,m+1) = sum(sum(temp.* covcell.ss{m+1}));
                end
                mlam = (lam'*pi/180)*(0:cols-1);
                f = cos(mlam)*cc*cos(mlam');
                f = f + sin(mlam)*ss*sin(mlam');
                % Saving the covariance fields in files
                if floor(theta(i)) < 10,
                    lat1 = ['00',num2str(theta(i),'%7.3f')];
                elseif floor(theta(i)) >= 10 && floor(theta(i)) < 100
                    lat1 = ['0',num2str(theta(i),'%7.3f')];
                elseif floor(theta(i)) >= 100
                    lat1 = num2str(theta(i),'%7.3f');
                end
                if i ~= j
                    if floor(theta(j)) < 10,
                        lat2 = ['00',num2str(theta(j),'%7.3f')];
                    elseif floor(theta(j)) >= 10 && floor(theta(i)) < 100
                        lat2 = ['0',num2str(theta(j),'%7.3f')];
                    elseif floor(theta(j)) >= 100
                        lat2 = num2str(theta(j),'%7.3f');
                    end
                elseif i == j
                    lat2 = lat1;
                    fvar(i,:) = diag(f)';
                end
                save([covpth, 'cov', lat1,'-', lat2, '.mat'], 'f')
                fprintf(fid,'%s \n',['cov', lat1,'-', lat2, '.mat']);
                cc = zeros(cols);
                ss = zeros(cols);
            end
            waitbar(i/nlat)
        end
        fclose(fid);
        fprintf('done\t')
    elseif ovartype == 1
        for i = 1:nlat
            for m = 0:lmax
                pnm = plm(m:lmax,m,trad(i));
                temp = pnm'*pnm;
                cc(m+1,m+1) = sum(sum(temp.* covcell.cc{m+1}));
                ss(m+1,m+1) = sum(sum(temp.* covcell.ss{m+1}));
            end
            mlam = (lam'*pi/180)*(0:cols-1);
            f = cos(mlam)*cc*cos(mlam');
            f = f + sin(mlam)*ss*sin(mlam');
            fvar(i,:) = diag(f)'; 
            cc = zeros(cols);
            ss = zeros(cols);
            waitbar(i/nlat)
        end
        fprintf('done\t')
    end
elseif ivartype == 2
	fprintf('diagonal of the covariance matrix ... ')
    covcell = covord(fname,pth,lmax,2);
    if delfile == 1, delete(fname), end

    if size(filcff,2) == 1
        for m = 1:lmax+1
            temp = filcff(m:end).^2;
            covcell.cc{m} = temp.*covcell.cc{m};
            covcell.ss{m} = temp.*covcell.ss{m};
        end
    elseif size(filcff,1)==size(filcff,2) && size(filcff,1)==lmax+1
        for m = 1:lmax+1
            temp = filcff(m:end,m).^2;
            covcell.cc{m} = temp.*covcell.cc{m};
            covcell.ss{m} = temp.*covcell.ss{m};
        end
    end
    if ovartype == 0
        fid = fopen([covpth,'cvflist.txt'],'W+');
        % Calculating the 2-D Fourier coefficients
        %
        for i = 1:nlat
            for j = i:nlat
                for m = 0:lmax
                    if i ~= j
                        pnmi = plm(m:lmax,m,trad([i;j]));
                        pnmj = pnmi(2,:);
                        pnmi(2,:) = [];
                    elseif i == j
                        pnmi = plm(m:lmax,m,trad(i));
                        pnmj = pnmi;
                    end
                    temp = pnmi.*pnmj;
                    cc(m+1,m+1) = sum(temp'.*covcell.cc{m+1});
                    ss(m+1,m+1) = sum(temp'.*covcell.ss{m+1});
                end
                mlam = (lam'*pi/180)*(0:cols-1);
                f = cos(mlam)*cc*cos(mlam');
                f = f + sin(mlam)*ss*sin(mlam');
                % Saving the covariance fields in files
                if floor(theta(i)) < 10,
                    lat1 = ['00',num2str(theta(i),'%7.3f')];
                elseif floor(theta(i)) >= 10 && floor(theta(i)) < 100
                    lat1 = ['0',num2str(theta(i),'%7.3f')];
                elseif floor(theta(i)) >= 100
                    lat1 = num2str(theta(i),'%7.3f');
                end
                if i ~= j
                    if floor(theta(j)) < 10,
                        lat2 = ['00',num2str(theta(j),'%7.3f')];
                    elseif floor(theta(j)) >= 10 && floor(theta(i)) < 100
                        lat2 = ['0',num2str(theta(j),'%7.3f')];
                    elseif floor(theta(j)) >= 100
                        lat2 = num2str(theta(j),'%7.3f');
                    end
                elseif i == j
                    lat2 = lat1;
                    fvar(i,:) = diag(f)';
                end
                save([covpth, 'cov', lat1,'-', lat2, '.mat'], 'f')
                fprintf(fid,'%s \n',['cov', lat1,'-', lat2, '.mat']);
                cc = zeros(cols);
                ss = zeros(cols);
            end
            waitbar(i/nlat)
        end
        fclose(fid);
        fprintf('done\t')
    elseif ovartype == 1
        for i = 1:nlat
            for m = 0:lmax

                pnm = plm(m:lmax,m,trad(i));
                temp = (pnm.^2)';
                cc(m+1,m+1) = sum(temp.*covcell.cc{m+1});
                ss(m+1,m+1) = sum(temp.*covcell.ss{m+1});
            end
            mlam = (lam'*pi/180)*(0:cols-1);
            f = cos(mlam)*cc*cos(mlam');
            f = f + sin(mlam)*ss*sin(mlam');
            fvar(i,:) = diag(f)'; 
            cc = zeros(cols);
            ss = zeros(cols);
            waitbar(i/nlat)
        end
        fprintf('done\t')
    end
end
close(hwb);

if ovartype == 0
    varargout{1} = 'cvflist.txt';
end

thetaRAD = theta * pi/180;
lamRAD 	 = lam * pi/180;

toc
