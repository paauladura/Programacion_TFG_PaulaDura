function [T,Tr,Trr,Trrr,Trrrr] = pshs(field,lambdaRAD,phiRAD,r,const,wbh,jflag)

% function determine the radial derivatives for the disturbing potential 
% T along the orbit. Also this file use spherical harmonics for
% calculations it is especially designed for use with series of latitude
% and longitude. So not a field is created but the derivates e.g. along a
% orbit.
%
% IN:
%    field ....... a priori gravity field in cs or sc format         [n,m]  
%    lambdaRAD ... longitude in [rad]                                [n,1]   
%    phiRAD ...... latitude in [rad]                                 [n,1]   
%    r ........... radius in [m]                                     [n,1]   
%    const ....... constant factors for GM and ae                    [2,1]  
%    wbh ......... waitbar handle                                    [1,1]  
%    jflag ....... reduce normal field (default: true)               [1,1]  
%
% OUT:
%    T ........... distrubing potentia                               [n,1]   
%    Tr .......... first radial derivative                           [n,1]  
%    Trr ......... second radial derivative                          [n,1]  
%    Trrr ........ third radial derivative                           [n,1]  
%    Trrrr ....... fourth radial derivative                          [n,1]  
%
% USES:
  %    cs2sc, normalklm, Legendre_mex, uberall/constants, uberall/checkcoor
%
% REMARKS:
%    r can be scalar

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
%   2021-04-12: MA, remove revision statement on deprecated function 
%   2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users)
%   2013-02-13: MR, change function names, brush up comments
%   2013-01-30: MA, comments
%   2013-01-23: MA, input in radian
%   2004-12-12: MW, change order of output arguments
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

%% INPUT CHECK and PREPARATION
if nargin < 7 || isempty(jflag), jflag  = true; end
if nargin < 6 || isempty(wbh),   wbh    = [];   end
if nargin < 5 || isempty(const), const  = [];   end

legmex = exist('Legendre_mex', 'file') == 3; % check if compiled Legendre_mex exists

% load necessary constants
if isempty(const)
    constants;
elseif numel(const) == 2
    GM = const(1);
    ae = const(2);
else
    error('CONST must be a 2x1 vector or empty.')
end

% Size determination for field
[row, col] = size(field);
if (row ~= col) && (col ~= 2*row-1), error('Input ''field'' not in cs or sc format'); end
field = cs2sc(field);
lmax  = row - 1;

% prepare the coordinates
[lambdaRAD, phiRAD, r] = checkcoor(lambdaRAD, phiRAD, r, ae, 'pointwise');
theta = (pi/2 - phiRAD);
didx   = 1:numel(lambdaRAD);

%% PREPARATION
T = NaN(size(lambdaRAD));
if nargout > 1, Tr    = NaN(size(lambdaRAD)); end
if nargout > 2, Trr   = NaN(size(lambdaRAD)); end
if nargout > 3, Trrr  = NaN(size(lambdaRAD)); end
if nargout > 4, Trrrr = NaN(size(lambdaRAD)); end

% substract reference field
if jflag, field = field - cs2sc(normalklm(lmax,'wgs84')); 
    warning('A reference field (WGS84) is removed from your coefficients')
end

% cutoff data in parts of 5000 elements
maxlength = 5000;

% prepare index
idx   = all(~isnan([lambdaRAD theta r]),2);
didx  = didx(idx);
lambdaRAD   = lambdaRAD(idx);
theta   = theta(idx);
r     = r(idx);

% calculation
if ~isempty(wbh), waitbar(0,wbh); end
parts = ceil(numel(lambdaRAD)/maxlength);
for I = 1:parts
%     if I == parts
%         idx = (I-1)*maxlength+1:numel(lambdaRAD);
%     else
%         idx = (I-1)*maxlength+1:I*maxlength;
%     end
    N = min([I*maxlength; numel(lambdaRAD)]);
    idx = (I-1)*maxlength+1:N;
    
    % prepare cosine and sine --> cos(m*lam) and sin(m*lam)
    l       = 0:row-1;
    m       = l';
    mlam    = m*lambdaRAD(idx)';         % matrix of size(length(m),length(lam))
    cosinus = cos(mlam);
    sinus   = sin(mlam);
    
    % prepare ratio Earth radius over radius
    n  = m(:,ones(length(idx),1))';
    TF = ae./r(idx) * ones(1,row);      % matrix of size(length(r),length(n))
    
    % prepare factors
    TK     = GM./r(idx);
    TFn    = TF.^n;
    if nargout >= 2, TrK    = -GM./r(idx).^2; TrF    = (n+1).*(TF.^n); end      % factors for Tr 
    if nargout >= 3, TrrK   = -TrK./r(idx);   TrrF   = (n+2).*TrF;     end      % factors for Trr 
    if nargout >= 4, TrrrK  = -TrrK./r(idx);  TrrrF  = (n+3).*TrrF;    end      % factors for Trrr 
    if nargout == 5, TrrrrK = -TrrrK./r(idx); TrrrrF = (n+4).*TrrrF;   end      % factors for Trrrr 
    
    %% CALCULATION
    TA = zeros(numel(idx),row); TB = zeros(numel(idx),row);
    if nargout >=2, TrA    = zeros(numel(idx),row); TrB    = zeros(numel(idx),row); end
    if nargout >=3, TrrA   = zeros(numel(idx),row); TrrB   = zeros(numel(idx),row); end
    if nargout >=4, TrrrA  = zeros(numel(idx),row); TrrrB  = zeros(numel(idx),row); end
    if nargout >=5, TrrrrA = zeros(numel(idx),row); TrrrrB = zeros(numel(idx),row); end
    for m = 0:row-1
        if m==0
            Cnm = field(:,row+m);               % get column with order 0
            if legmex                           % fully normalized Legendre Polynoms
                P = Legendre_mex(l, m, theta(idx));
            else
                P = plm(l, m, theta(idx));  
            end
            % ------------------------------------------------------------------
            TP      = TFn.*P;             % calc for fourth radial derivative 
            TA(:,1) = TP*Cnm;
            % ----------------------------
            if nargout >= 2
                TrP         = TrF.*P;     % calc for first radial derivative 
                TrA(:,1)    = TrP*Cnm;
            end
            % ----------------------------
            if nargout >= 3
                TrrP        = TrrF.*P;    % calc for second radial derivative 
                TrrA(:,1)   = TrrP*Cnm;
            end
            % -----------------------------
            if nargout >= 4
                TrrrP       = TrrrF.*P;   % calc for third radial derivative 
                TrrrA(:,1)  = TrrrP*Cnm;
            end
            % ------------------------------------------------------------------
            if nargout == 5
                TrrrrP      = TrrrrF.*P;  % calc for fourth radial derivative 
                TrrrrA(:,1) = TrrrrP*Cnm;
            end
            % ------------------------------------------------------------------
        else
            Cnm = field(m+1:end,row+m);        % get Cnm coefficients for order m
            Snm = field(m+1:end,row-m);        % get Snm coefficients for order m
            if legmex                           % fully normalized Legendre Polynoms
                P = Legendre_mex(l, m, theta(idx));
            else
                P = plm(l, m, theta(idx));  
            end
            P = P(:,m+1:end);

            % ------------------------------------------------------------------
            TP     = TFn(:,m+1:end).*P;        % calc for potential
            TA(:,m+1) = TP*Cnm;
            TB(:,m+1) = TP*Snm;
            % ------------------------------------------------------------------
            if nargout >= 2
                TrP           = TrF(:,m+1:end).*P;        % calc for first radial derivative 
                TrA(:,m+1)    = TrP*Cnm;
                TrB(:,m+1)    = TrP*Snm;
            end
            % ------------------------------------------------------------------
            if nargout >= 3
                TrrP          = TrrF(:,m+1:end).*P;   % calc for second radial derivative 
                TrrA(:,m+1)   = TrrP*Cnm;
                TrrB(:,m+1)   = TrrP*Snm;
            end
            % ------------------------------------------------------------------
            if nargout >= 4
                TrrrP         = TrrrF(:,m+1:end).*P;  % calc for third radial derivative 
                TrrrA(:,m+1)  = TrrrP*Cnm;
                TrrrB(:,m+1)  = TrrrP*Snm;
            end
            % ------------------------------------------------------------------
            if nargout == 5
                TrrrrP        = TrrrrF(:,m+1:end).*P; % calc for fourth radial derivative 
                TrrrrA(:,m+1) = TrrrrP*Cnm;
                TrrrrB(:,m+1) = TrrrrP*Snm;
            end
            % ------------------------------------------------------------------
            
        end
    end
    
    % now do the final sumation
    T(didx(idx)) = TK.*sum(TA.*cosinus'+TB.*sinus',2);
    if nargout >= 2, Tr(didx(idx))    = TrK.*sum(TrA.*cosinus'+TrB.*sinus',2);          end;
    if nargout >= 3, Trr(didx(idx))   = TrrK.*sum(TrrA.*cosinus'+TrrB.*sinus',2);       end;
    if nargout >= 4, Trrr(didx(idx))  = TrrrK.*sum(TrrrA.*cosinus'+TrrrB.*sinus',2);    end;
    if nargout == 5, Trrrr(didx(idx)) = TrrrrK.*sum(TrrrrA.*cosinus'+TrrrrB.*sinus',2); end;
    
    % CleanUp
    clear TrP TrF TrA TrB TrK TrrP TrrF TrrA TrrB TrrK TrrrP TrrrF TrrrA TrrrB TrrrK;
    clear TrrrrP TrrrrF TrrrrA TrrrrB TrrrrK TP TFn TA TB TK P Cnm Snm cosinus sinus;
    clear l m n mlam idx;
    
    if ~isempty(wbh), waitbar(I/parts,wbh); end
    
end % end for I=1:parts

