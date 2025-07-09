function f = gshs_grid(field, lamRAD, phiRAD, a_E, varargin)

% GSHSAG calculates a global spherical harmonic synthesis for any grid 
% defined by lam and phi (both vectors). The radius must be scalar.
% Multiple sets of Stokes coefficients can be inserted as 3D matrix. 
% 
% f = gshs_grid(field, lamRAD, phiRAD, a_E)
%
% IN:
%    field ....... gravity field in |c\s| or /s|c\ format            [r, c, d]
%                  (in case of multiple SH coefficients: one set per layer d of the 3D matrix)
%    lamRAD ...... longitude [rad]                                   [n, 1]
%    phiRAD ...... latitude  [rad]                                   [m, 1]
%    a_E ......... semi major axis of the Earth rotational ellipsoid [1, 1]
% OPTIONAL:
%    'height' .... (default: 0), height [m]                        [scalar]
%    'max_lm' .... maximum degree/order (default: determined from field)
%                                                                  [scalar]
%    'quant' ..... optional argument, defining the field quantity: [string]
%                  - 'potential' ... (default), potential [m^2/s^2], needs 'GM'
%                  - 'tr' .......... gravity disturbance [mGal], needs 'GM'
%                  - 'trr' ......... 2nd rad. derivative [E], needs 'GM'
%                  - 'none' ........ coefficients define the output
%                  - 'geoid' ....... geoid height [m] 
%                  - 'dg', 'gravity' gravity anomaly [mGal], needs 'GM'
%                  - 'slope' ....... slope [arcsec]
%                  - 'water' ....... equivalent water thickness [m] 
%                  - 'smd' ......... surface mass density
%    'sub_WGS84' . (default: true) determines whether WGS84 is subtracted.
%    'GM' ........ geocentric gravitational constant GM
%    'legendre' .. (default: 'plm') Legendre functions algorithm:
%                  - 'plm' ... unstable (for d/o > ~1800) PLM
%                  - 'mex' ... fast, X-number stabilized LEGENDRE_MEX
%    'waitbar' ... if set, shows a waitbar (default: false)
%    'curvature' . (default: false) consider curvature of the chosen quantity [bool]
%
% OUT:
%    f ....... field quantity                                       [n, m]
%
% EXAMPLE: see SHBUNDLE/example/example_gshs_grid.m
%          see SHBUNDLE/example/example_slepian.m
%
% USES:
%    vec2cs, normalklm, eigengrav, plm, Legendre\_mex, cs2sc, sc2cs, 
%    UBERALL/checkcoor, UBERALL/twaitbar
%
% SEE ALSO:
%    GSHS_, GSHS_PTW, REGIONALSHS

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), DoGE, UofC
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2021-10-20: MA, possibility of multiple sets of SH coefficients
%    2021-04-12: MA, remove comment on deprecated function 
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users; avoidable via getopt.m)
%    2015-05-22: MR, reprogram parameter interface (hence, rename function), 
%                    revise help text and code
%    2014-01-15: MR, revise help text, beautify code
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments/removing of smoothing option
%    2013-01-23: MA, input of plm in radian
%    2012-03-06: MW, add curvature calculations
%    2009-03-04: MW, change input co-latitude -> latitude
%    2005-10-25: MW, initial version
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

%--------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%--------------------------------------------------------------------------


%% define defaults and check optional parameters
defaultparams = {'height', 0;
                 'max_lm', inf;
                 'quant', 'potential';
                 'sub_wgs84', true;
                 'gm', 1;
                 'legendre', 'plm';
                 'waitbar', false;
                 'curvature', false}; 
params = getopt(defaultparams, false, varargin);  

% check Legendre function algorithm
plm_func = false;
switch params.legendre
    case 'mex'
        if exist('Legendre_mex', 'file') == 3; % check if compiled Legendre_mex exists
            plm_func = @Legendre_mex;
        else
            warning('Legendre_mex is not compiled! Falling back to plm...');
            plm_func = @plm;
        end
    case 'plm'
        plm_func = @plm;
    otherwise
        error('Unknown Legendre function.');
end

% field size determination, rearrange field and subtract reference field
[row, col,depth] = size(field);

if col == row              % field in |C\S|-format 
    for dd = depth:-1:1
        copy(:,:,dd) = cs2sc(field(:,:,dd));  % convert to /S|C\-format
    end
    field = copy;
    col = 2* row -1;
elseif col ~= 2 * row - 1
   error('Input "field" not in cs or sc format');
% else: field is already in /S|C\-format
end

if params.max_lm > (row - 1)     % desired max_lm is bigger than what the field provides
    params.max_lm = row - 1;     % adjust max_lm 
elseif params.max_lm < (row - 1) % if max_lm is smaller than what the field provides
    field = field(1:(params.max_lm + 1), (row - params.max_lm):(row + params.max_lm),:); % adjust field
    row = size(field, 1); % update row
    col = 2* row -1;
% else: everything is ok
end

if params.sub_wgs84
    field = bsxfun(@minus, field, cs2sc(full(normalklm(params.max_lm, 'wgs84'))));
    if params.display_warning == true
        warning('A reference field (WGS84) is removed from your coefficients')
    end
end

% check if r, lamRAD and phiRAD are vectors
[lamRAD, phiRAD, ~] = checkcoor(lamRAD, phiRAD, 0, 0, 'grid');
theRAD              = (pi/2 - phiRAD);
if ~isscalar(params.height)
    error('''height'' must be scalar.');
end

% prepare l, m, often used variables
m       = (0:params.max_lm);
l       = m';
mlam    = (lamRAD * m)';           
len_the = length(theRAD);

% apply transfer function
transf  =  eigengrav(l, params.quant, params.height, [params.gm, a_E]);
field   =  reshape(diag(transf(:)) * field, row, col, depth);
 
%----------------------------------------------------------------------------
% CALCULATION
%----------------------------------------------------------------------------
if params.waitbar
    WAIT = twaitbar('init', [], 'gshs_grid: calculating');
end;

if params.curvature
    TAp  = []; TAp(len_the, row,depth) = 0; % faster than TAp = zeros(len_the, row)
    TBp  = []; TBp(len_the, row,depth) = 0;
    TAl  = []; TAl(len_the, row,depth) = 0;
    TBl  = []; TBl(len_the, row,depth) = 0;
    TApp = []; TApp(len_the, row,depth) = 0; 
    TBpp = []; TBpp(len_the, row,depth) = 0;
    TAll = []; TAll(len_the, row,depth) = 0; 
    TBll = []; TBll(len_the, row,depth) = 0;
    TApl = []; TApl(len_the, row,depth) = 0; 
    TBpl = []; TBpl(len_the, row,depth) = 0;
    
    % unwrapping to avoid if ... else for order m = 0 (i.e. Snm = 0)
    Cnm = squeeze(field(:, row,:));              % get Cnm coefficients for order 0
    [P, dP, ddP] = plm_func(l, 0, theRAD); % calc fully normalized Legendre Polynoms
    TAp(:, 1,:)  = -dP * Cnm; % all variables multiplied with Snm are 0, of course  
    TBl(:, 1,:)  =  -P * Cnm;
    TApp(:, 1,:) = ddP * Cnm;
    TAll(:, 1,:) =  -P * Cnm;
    TBpl(:, 1,:) =  dP * Cnm;
    for m = 1:(row - 1)
        Cnm = squeeze(field(:, row + m,:));          % get Cnm coefficients for order m
        Snm = squeeze(field(:, row - m,:));          % get Snm coefficients for order m
        [P, dP, ddP] = plm_func(l, m, theRAD); % calc fully normalized Legendre Polynoms
        TAp(:, m+1,:)  = -dP * Cnm;  TBp(:, m+1,:)  = -dP * Snm;
        TAl(:, m+1,:)  =   P * Snm;  TBl(:, m+1,:)  =  -P * Cnm;
        TApp(:, m+1,:) = ddP * Cnm;  TBpp(:, m+1,:) = ddP * Snm;
        TAll(:, m+1,:) =  -P * Cnm;  TBll(:, m+1,:) =  -P * Snm;
        TApl(:, m+1,:) = -dP * Snm;  TBpl(:, m+1,:) =  dP * Cnm;
        if params.waitbar; WAIT = twaitbar(m / (row - 1), WAIT); end;
    end
    
    % sum the partial derivatives

    cosmlam = cos(mlam);
    sinmlam = sin(mlam);
    for dd = depth:-1:1 %
        fp  = TAp(:,:,dd)  * cosmlam + TBp(:,:,dd)  * sinmlam; fp2 = fp.^2;
        fl  = TAl(:,:,dd)  * cosmlam + TBl(:,:,dd)  * sinmlam; fl2 = fl.^2;
        fpp = TApp(:,:,dd) * cosmlam + TBpp(:,:,dd) * sinmlam;
        fll = TAll(:,:,dd) * cosmlam + TBll(:,:,dd) * sinmlam;
        fpl = TApl(:,:,dd) * cosmlam + TBpl(:,:,dd) * sinmlam;
    
        % now do the final summation
        f(:,:,dd) = ((1 + fp2) .* fll - 2 .* fp .* fl .* fpl + (1 + fl2) .* fpp) ./ (2 .* (realsqrt(1 + fp2 + fl2)).^3);
    end
else
    TA = []; TA(len_the, row,depth) = 0; % faster than TA = zeros(length(idx), row);
    TB = []; TB(len_the, row,depth) = 0;   

    % unwrapping to avoid if ... else for order m = 0
    Cnm = squeeze(field(:, row,:));                % get Cnm coefficients for order 0
    TA(:, 1,:) = plm_func(l, 0, theRAD) * Cnm; % for m = 0: Snm = 0, hence TB also = 0
    for m = 1:row-1
        Cnm = squeeze(field(:, row + m,:));     % get Cnm coefficients for order m
        Snm = squeeze(field(:, row - m,:));     % get Snm coefficients for order m
        P = plm_func(l, m, theRAD);  % calc fully normalized Legendre Polynoms
        % ------------------------------------------------------------------
        TA(:, m + 1,:) = P * Cnm;
        TB(:, m + 1,:) = P * Snm;
        if params.waitbar; WAIT = twaitbar(m / (row - 1), WAIT); end;
    end
    % now do the final summation
    cosmlam = cos(mlam);
    sinmlam = sin(mlam);
    for dd = depth:-1:1 
        f(:,:,dd) = TA(:,:,dd) * cosmlam  + TB(:,:,dd) * sinmlam;
    end
end

if params.waitbar; twaitbar('close', WAIT); end;

