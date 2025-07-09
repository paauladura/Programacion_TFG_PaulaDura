function [tensV,gradV] = tenspshs(field,lamRAD,phiRAD,r,lmax,const,drtype,wbh,jflag)

% TENSPSHS determines the tensor and gradient of the potential pointwise at
% the location (r,\phi,\lambda)
%
% IN:
%    field ..... a priori gravity field in cs or sc format                  [n x m]  
%    lamRAD .... longitude in [rad]                                         [t x 1]  
%    phiRAD .... latitude in [rad]                                          [t x 1]  
%    r ......... radius in [m]                                              [t x 1]  
%    lmax ...... maximum degree (default: inf)                              [1 x 1]  
%    const ..... constants GM (product of gravitational constant G          [2 x 1]  
%                with mass of the central body M) and ae (semi-major
%                axis of the central body)
%    drtype .... defines the type of output                                 [string]  
%                - 'deriv' are the derivatives towards r, \phi (p), \lambda (l)
%                - 'xyz'   are the derivatives towards {x,y,z}
%                          but in the cartesian frame
%                - 'lscs'  (default) are the derivatives towards {x,y,z}
%                          but in the local spherical coordinate frame
%                - [t x 9] rotation matrix to rotate the results into an 
%                          arbitrary coordinate frame
%    wbh ....... waitbar handle                                             [1 x 1]  
%    jflag ..... 1 = (default) reference field is subtracted                [1 x 1]  
%                0 = reference field is not subtracted
% OUT: 
%    tensV ..... tensor in [1/s^2]:                                         [t x 9]  
%                            - 'deriv':
%                              [ Vrr Vrp Vrl Vpr Vpp Vpl Vlr Vlp Vll ]
%                            - 'xyz', 'lscs': 
%                              [ Vxx Vxy Vxz Vyx Vyy Vyz Vxz Vyz Vzz ]
%    gradV ..... gradient in [m/s^2]             [t x 3]  
%                            - 'deriv':
%                              [ Vr Vp Vl ]
%                            - 'xyz', 'lscs': 
%                              [ Vx Vy Vz ]
%
% USES:
%    normalklm, plm, cs2sc, sc2cs, uberall/multmat
%    uberall/multmatvek, uberall/constants, uberall/checkcoor

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), GI, Uni Stuttgart
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users)
%    2013-02-26: MR, help text: adding description of functionality, code:
%                    adding input check for variable 'drtype'
%    2013-02-25: MR, clarification of meanings in help text
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments
%    2013-01-23: MA, input in radian
%    2010-07-16: MW, initial version
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


%% INPUT CHECK
if nargin < 9 || isempty(jflag),  jflag  = 1;      end
if nargin < 8 || isempty(wbh),    wbh    = [];     end
if nargin < 7 || isempty(drtype), drtype = 'lscs'; end
if nargin < 6 || isempty(const),  const  = [];     end

if ~exist('lmax', 'var') || isempty(lmax)
    lmax = inf;
end

% check parameter in drtype
if ~strcmpi(drtype,'deriv') && ~strcmpi(drtype,'xyz') && ~strcmpi(drtype,'lscs') && ~any(size(drtype) == 9)
    error('unknown value for DRTYPE');
end
   
% load necessary constants
if isempty(const)
    constants;
elseif numel(const) == 2
    GM = const(1);
    ae = const(2);
else
    error('CONST must be a 2x1 vector or empty.')
end

% check the rotation matrix

%% PREPARATION
% prepare the coordinates
[lamRAD,phiRAD,r] = checkcoor(lamRAD,phiRAD,r,ae,'pointwise');
if ~ischar(drtype) && any(size(drtype) ~= [numel(lamRAD) 9])
    error('Rotation matrix must have the same number of rows as LAM and 9 columns');
end

% preparing the input field
[row col] = size(field);
if (row~=col) && (col~=2*row-1)
   error('Input FIELD not in cs or sc format');
elseif col==2*row-1
   field = sc2cs(field);
end
if (lmax == inf) || (row - 1 < lmax), lmax = row-1; end
field = field(1:lmax+1,1:lmax+1); 
field = cs2sc(field,0);
row   = size(field,1);

% substract reference field if requested
if jflag, field = field - cs2sc(normalklm(lmax,'wgs84')); 
    warning('A reference field (WGS84) is removed from your coefficients')
end

% prepare the output
tensV = NaN.*ones(length(lamRAD),9);
gradV = NaN.*ones(length(lamRAD),3);
didx  = 1:numel(lamRAD);
didx  = didx(all(~isnan([lamRAD phiRAD r]),2));

%% CALCULATION
maxlength = 512;
if ~isempty(wbh), waitbar(0,wbh); end

% process data (piecewise)
parts = ceil(numel(didx)/maxlength);
for I = 1:parts
    % select the data
    if I == parts
        idx = didx((I-1)*maxlength+1:end);
    else
        idx = didx((I-1)*maxlength+1:I*maxlength);
    end
    
    % prepare cosine and sine --> cos(m*lam) and sin(m*lam) and  get
    % co-latitude 
    theRAD      = (pi/2 - phiRAD(idx)); 
    m       = (0:row-1)';
    mlam    = m*lamRAD(idx)';             % matrix of size(length(m),length(lam))
    cosinus = cos(mlam);
    sinus   = sin(mlam);
    
    % prepare ratio Earth radius over radius
    l   =    m(:,ones(length(idx),1))';
    TF  =    ae./r(idx) * ones(1,row);     % matrix of size(length(r),length(n))
    TF  =    TF.^l;
    TFn1 = - TF.*(l+1);
    TFn2 =   TF.*(l+1).*(l+2);
    Tr   =   GM./r(idx);
    Trr  =   Tr./r(idx);                  
    Trrr =   Trr./r(idx);
    
    % m-loop: summation over the order
    TrA  = zeros(numel(idx),row);
    TrB  = zeros(numel(idx),row);
    TpA  = zeros(numel(idx),row);
    TpB  = zeros(numel(idx),row);
    TlA  = zeros(numel(idx),row);
    TlB  = zeros(numel(idx),row);
    TrrA = zeros(numel(idx),row);
    TrrB = zeros(numel(idx),row);
    TrpA = zeros(numel(idx),row);
    TrpB = zeros(numel(idx),row);
    TrlA = zeros(numel(idx),row);
    TrlB = zeros(numel(idx),row);
    TppA = zeros(numel(idx),row);
    TppB = zeros(numel(idx),row);
    TplA = zeros(numel(idx),row);
    TplB = zeros(numel(idx),row);
    TllA = zeros(numel(idx),row);
    TllB = zeros(numel(idx),row);
    for m = 0:row-1
        if m==0
            Cnm         = field(:,row+m);         % get column with order 0
            [P,dP,ddP]  = plm(m:lmax,m,theRAD);       % calc fully normalzied Legendre Polynoms 
            
            % Gradient ----------------------------------------------------
            TrP         = TFn1.*P;                % calc the radial derivative
            TrA(:,1)    = TrP*Cnm;
            TrB(:,1)    = zeros(size(TrA(:,1))); 
            % -----------------
            TpP         = TF.*-dP;                 % calc the derivative towards phi
            TpA(:,1)    = TpP*Cnm;
            TpB(:,1)    = zeros(size(TpA(:,1))); 
            % -----------------
            TlA(:,1)    = zeros(size(TpA(:,1)));  % the derivative towards lambda is zero for order zero
            TlB(:,1)    = zeros(size(TpA(:,1)));
            
            % Tensor ------------------------------------------------------
            TrrP        = TFn2.*P;                % calc second radial derivative
            TrrA(:,1)   = TrrP*Cnm;
            TrrB(:,1)   = zeros(size(TrrA(:,1)));  
            % ---------------------
            TrpP        = TFn1.*-dP;               % calc radial/derivative towards phi
            TrpA(:,1)   = TrpP*Cnm;
            TrpB(:,1)   = zeros(size(TrpA(:,1)));           
            % ---------------------
            TrlA(:,1)   = zeros(size(TrpA(:,1))); % the derivative towards lambda is zero for order zero
            TrlB(:,1)   = zeros(size(TrpA(:,1)));         
            % ---------------------
            TppP        = TF.*ddP;                % calc second derivative towards phi
            TppA(:,1)   = TppP*Cnm;
            TppB(:,1)   = zeros(size(TppA(:,1)));   
            % ---------------------
            TplA(:,1)   = zeros(size(TppA(:,1))); % the derivative towards lambda is zero for order zero
            TplB(:,1)   = zeros(size(TppA(:,1)));
            % ---------------------
            TllA(:,1)   = zeros(size(TppA(:,1))); % the derivative towards lambda is zero for order zero
            TllB(:,1)   = zeros(size(TppA(:,1)));

        else
            Cnm         = field(m+1:end,row+m);   % get Cnm coefficients for order m
            Snm         = field(m+1:end,row-m);   % get Snm coefficients for order m
            [P,dP,ddP]  = plm(m:lmax,m,theRAD);       % calc fully normalzied Legendre Polynoms 
            
            % Gradient ----------------------------------------------------
            TrP         = TFn1(:,m+1:end).*P;     % calc the radial derivative
            TrA(:,m+1)  = TrP*Cnm;
            TrB(:,m+1)  = TrP*Snm;
            % --------------------
            TpP         = TF(:,m+1:end).*-dP;      % calc the derivative towards phi
            TpA(:,m+1)  = TpP*Cnm;
            TpB(:,m+1)  = TpP*Snm;
            % --------------------
            TlP         = TF(:,m+1:end).*P;       % calc the derivative towards lambda 
            TlA(:,m+1)  = -m.*(TlP*Cnm);
            TlB(:,m+1)  =  m.*(TlP*Snm);

            % Tensor ------------------------------------------------------
            TrrP        = TFn2(:,m+1:end).*P;     % calc second radial derivative
            TrrA(:,m+1) = TrrP*Cnm;
            TrrB(:,m+1) = TrrP*Snm; 
            % ---------------------
            TrpP        = TFn1(:,m+1:end).*-dP;    % calc radial/derivative towards phi
            TrpA(:,m+1) = TrpP*Cnm;
            TrpB(:,m+1) = TrpP*Snm;            
            % ---------------------
            TrlP        = TFn1(:,m+1:end).*P;     % calc radial/derivative towards lambda
            TrlA(:,m+1) = -m.*(TrlP*Cnm);
            TrlB(:,m+1) =  m.*(TrlP*Snm);
            % ---------------------
            TppP        = TF(:,m+1:end).*ddP;     % calc second derivative towards phi
            TppA(:,m+1) = TppP*Cnm;
            TppB(:,m+1) = TppP*Snm;
            % ---------------------
            TplP        = TF(:,m+1:end).*-dP;      % calc derivative towards phi/lambda
            TplA(:,m+1) = -m.*(TplP*Cnm);
            TplB(:,m+1) =  m.*(TplP*Snm);
            % ---------------------
            TllP        = TF(:,m+1:end).*P;       % calc derivative towards phi/lambda
            TllA(:,m+1) = -m.^2.*(TllP*Cnm);
            TllB(:,m+1) = -m.^2.*(TllP*Snm);

        end
    end
    
    % Now do the final sumation
    % Gradient ----------------------------------------------------
    gradV(idx,1) = Trr .*sum(TrA .*cosinus' + TrB.*sinus'  , 2);
    gradV(idx,2) = Tr  .*sum(TpA .*cosinus' + TpB.*sinus'  , 2);
    gradV(idx,3) = Tr  .*sum(TlA .*sinus'   + TlB.*cosinus', 2);
    % Tensor ------------------------------------------------------
    tensV(idx,1) = Trrr.*sum(TrrA.*cosinus' + TrrB.*sinus',  2);
    tensV(idx,2) = Trr .*sum(TrpA.*cosinus' + TrpB.*sinus',  2);
    tensV(idx,3) = Trr .*sum(TrlA.*sinus'   + TrlB.*cosinus',2);
    tensV(idx,4) = tensV(idx,2);
    tensV(idx,5) = Tr  .*sum(TppA.*cosinus' + TppB.*sinus',  2);
    tensV(idx,6) = Tr  .*sum(TplA.*sinus'   + TplB.*cosinus',2);    
    tensV(idx,7) = tensV(idx,3);
    tensV(idx,8) = tensV(idx,6);
    tensV(idx,9) = Tr  .*sum(TllA.*cosinus' + TllB.*sinus',  2);
    
    % prepare output and rotate if desired
    if strcmpi(drtype,'xyz') || strcmpi(drtype,'lscs') || any(size(drtype) == 9)
        % Tensor ------------------------------------------------------
        tensV(idx,2) = tensV(idx,2)./r(idx) - gradV(idx,2)./(r(idx).^2);
        tensV(idx,3) = tensV(idx,3)./r(idx)./cos(phiRAD(idx)) - gradV(idx,3)./(r(idx).^2)./cos(phiRAD(idx));
        tensV(idx,4) = tensV(idx,2);
        tensV(idx,5) = tensV(idx,5)./(r(idx).^2) + gradV(idx,1)./r(idx);
        tensV(idx,6) = tensV(idx,6)./(r(idx).^2)./cos(phiRAD(idx)) + gradV(idx,3).*tan(phiRAD(idx))./(r(idx).^2)./cos(phiRAD(idx));
        tensV(idx,7) = tensV(idx,3);
        tensV(idx,8) = tensV(idx,6);
        tensV(idx,9) = tensV(idx,9)./(r(idx).^2)./(cos(phiRAD(idx)).^2) + gradV(idx,1)./r(idx) - gradV(idx,2)./(r(idx).^2).*tan(phiRAD(idx));
        % Gradient ----------------------------------------------------
        gradV(idx,2) = gradV(idx,2)./r(idx);
        gradV(idx,3) = gradV(idx,3)./(r(idx).*cos(phiRAD(idx)));
    end
    
    if strcmpi(drtype,'xyz') || any(size(drtype) == 9)
        cp = cos(phiRAD(idx));
        sp = sin(phiRAD(idx));
        cl = cos(lamRAD(idx));
        sl = sin(lamRAD(idx));
        J  = [cp.*cl, -sp.*cl, -sl, cp.*sl, -sp.*sl, cl, sp, cp, zeros(size(sp))];
        gradV(idx,:) = multmatvek(J,gradV(idx,:));
        tensV(idx,:) = multmat(J,multmat(tensV(idx,:),J(:,[1 4 7 2 5 8 3 6 9])));
    end
    
    if any(size(drtype)==9)
        gradV(idx,:) = multmatvek(drtype(idx,:),gradV(idx,:));
        tensV(idx,:) = multmat(drtype(idx,:),multmat(tensV(idx,:),drtype(idx,[1 4 7 2 5 8 3 6 9])));
    end
    
    % CleanUp
    clear Tr Trr Trrr TF TFn1 TFn2 TrA TrB TrP TpA TpB TpP TlA TlB TlP;
    clear TrrA TrrB TrrP TrpA TrpB TrpP TrlA TrlB TrlP TppA TppB TppP;
    clear TplA TplB TplP TllA TllB TllP P dP ddP;
    clear m l mlam theRAD idx Cnm Snm cosinus sinus;
    
    if ~isempty(wbh), waitbar(I/parts,wbh); end
end % end for I=1:parts


