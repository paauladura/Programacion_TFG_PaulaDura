function [Tr,Tth,Tlam,T] = nablaPot(field, lambdaRAD, thetaRAD, r)

% function determine the radial derivatives for the disturbing potential 
% T along the orbit. Also this file use spherical harmonics for
% calculations it is especially designed for use with series of latitude
% and longitude. So not a field is created but the derivates e.g. along a
% orbit.
%
% IN:
%    field ......... a priori gravity field in cs or sc format     [n x m] 
%    lambdaRAD ..... longitude in [rad]                            [n x 1]   
%    thetaRAD ...... co-latitude in [rad]                          [n x 1]   
%    r ............. radius in [m]                                 [n x 1]  
%
% OUT:
%    T ............ distrubing potential                           [n x 1]  
%    Tr ........... first radial derivative                        [n x 1]  
%    Tth .......... first latitudal derivative                     [n x 1]  
%    Tlam ......... first longitudal derivative                    [n x 1]  
%
% EXAMPLE: see SHBUNDLE/examples/example_nablaPot.m
%
% USES:
%  vec2cs, normalklm, cs2sc, plm, uberall/multmatvek, uberall/checkcoor, 
%  uberall/constants
% 
% REMARKS:
%    r can be scalar

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%   Matthias WEIGELT (MW), DoGE, UofC  
%   Markus ANTONI (MA), GI, Uni Stuttgart
%   Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users)
%   2014-01-13: MR, remove unnecessary code, beautify code, update example
%   2013-02-13: MR, change function names, brush up comments
%   2013-01-29: MA, comments
%   2013-01-23: MA, input of plm in radian
%   2003-07-22: MW, initial version
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

% -------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
% -------------------------------------------------------------------------
% load necessary constants
constants;

% Size determination for field
[row, col] = size(field);
if (row ~= col) && (col ~= 2 * row - 1)
   error('Input ''field'' not in cs or sc format');
elseif col ~= 2 * row - 1
    % if data is in |C\S|-format we transfer it to /S|C\-format
   field = cs2sc(field);
   [row, ~] = size(field); %Hmoh ~ = col
end
lmax = row - 1;

% prepare the coordinates
[lambdaRAD, phiRAD, r] = checkcoor(lambdaRAD, pi/2 - thetaRAD, r, ae, 'pointwise');
thetaRAD = pi/2 - phiRAD;

%----------------------------------------------------------------------------
% PREPARATION
%----------------------------------------------------------------------------
T        = NaN.*ones(size(lambdaRAD));
Tr       = NaN.*ones(size(lambdaRAD));
Tth      = NaN.*ones(size(lambdaRAD));
Tlam     = NaN.*ones(size(lambdaRAD));

[idx, ~] = find(~isnan([lambdaRAD thetaRAD r])); % returns row/col indices (not linear), but col indices are not used
idx      = unique(idx);

lambdaRAD   = lambdaRAD(idx);
thetaRAD   = thetaRAD(idx);
r        = r(idx);

% substract reference field
field = field - cs2sc(normalklm(lmax, 'wgs84'));
warning('A reference field (WGS84) is removed from your coefficients')

% cutoff data in parts of 1500 elements
maxlength = 1500;

if length(lambdaRAD) > maxlength
    parts = ceil(length(lambdaRAD) / maxlength);
    hidx  = idx;
    hlam  = lambdaRAD;
    hth   = thetaRAD;
    hr    = r;
    jflag = 1;
else
    parts = 1;
    jflag = 0;
end
h = waitbar(0, 'Please wait...');

% calculation 
for I = 1:parts
    waitbar(I / parts, h)
    
    % distribute parts of input vectors
    if jflag == 1  
        if I == parts
            idx    = hidx(((I-1) * maxlength + 1):end);
            lambdaRAD = hlam(((I-1) * maxlength + 1):end);
            thetaRAD = hth(((I-1) * maxlength + 1):end);
            r      = hr(((I-1) * maxlength + 1):end);
        else
            idx    = hidx(((I-1) * maxlength + 1):(I * maxlength));
            lambdaRAD = hlam(((I-1) * maxlength + 1):(I * maxlength));
            thetaRAD = hth(((I-1) * maxlength + 1):(I * maxlength));
            r      = hr(((I-1) * maxlength + 1):(I * maxlength));
        end
    end
    idxlen = length(idx);
    
    % prepare cosine and sine --> cos(m*lam) and sin(m*lam)
    m       = (0:row - 1)';
    mlam    = m * lambdaRAD';          % matrix of size(length(m), length(lam))
    cosinus = cos(mlam);
    sinus   = sin(mlam);
    dcos    = -m(:, ones(1, length(lambdaRAD))) .* sinus;
    dsin    =  m(:, ones(1, length(lambdaRAD))) .* cosinus;
    
    % prepare ratio Earth radius over radius
    l  = m';
    n  = m(:, ones(length(r), 1))';
    TF = ae ./ r * ones(1, row);    % matrix of size(length(r),length(n))
    
    % prepare factors for \nabla_r T
    TrK = -GM ./ r .^ 2;         % factors for Tr 
    TrF = (n + 1) .* (TF .^ n);
    
    % prepare factors for T, \nabla_\thetaRAD T and \nabla_\lambdaRAD T
    TK  = GM ./ r;               % factors for T
    TFn = TF .^ n;
    
    % initialize with zero
    TrA   = zeros(idxlen, row);
    TrB   = zeros(idxlen, row);
    TthA  = zeros(idxlen, row);
    TthB  = zeros(idxlen, row);
    TlamA = zeros(idxlen, row);
    TlamB = zeros(idxlen, row);
    TA    = zeros(idxlen, row);
    TB    = zeros(idxlen, row);

    %----------------------------------------------------------------------------
    % CALCULATION
    %----------------------------------------------------------------------------
    for m = 0:row-1
        Cnm = field(:, row + m);            % get Cnm coefficients for order m
        if m == 0
            Snm = zeros(row, 1);            % there are no Sn0 coefficients
        else
            Snm = field(:, row - m);        % get Snm coefficients for order m
        end
        [P, dP] = plm(l, m, thetaRAD);      % calc fully normalzied Legendre Polynoms
        % ------------------------------------------------------------------
        TrP                  = TrF .* P;    % calc for first radial derivative 
        TrA(:, m+1)   = TrP * Cnm;
        TrB(:, m+1)   = TrP * Snm;
        % ------------------------------------------------------------------
        TthP                 = TFn .* -dP;   % calc for latitudal derivative
        TthA(:, m+1)  = TthP * Cnm;
        TthB(:, m+1)  = TthP * Snm;
        % ------------------------------------------------------------------
        TlamP                = TFn .* P;    % calc for longitudal derivative
        TlamA(:, m+1) = TlamP * Cnm;
        TlamB(:, m+1) = TlamP * Snm;    
        % ------------------------------------------------------------------
        TP                   = TFn .* P;    % calc for potential
        TA(1:idxlen, m+1)    = TP * Cnm;
        TB(1:idxlen, m+1)    = TP * Snm;
    end 
    
    % final sumation
    Tr(idx)    = TrK .* sum(TrA   .* cosinus' + TrB   .* sinus', 2);
    Tth(idx)   = TK  .* sum(TthA  .* cosinus' + TthB  .* sinus', 2);
    Tlam(idx)  = TK  .* sum(TlamA .* dcos'    + TlamB .* dsin',  2);
    T(idx)     = TK  .* sum(TA    .* cosinus' + TB    .* sinus', 2); 
end % for I = 1:parts 
close(h)

