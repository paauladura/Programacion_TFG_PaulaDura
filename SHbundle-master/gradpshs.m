function gradV = gradpshs(field,lamRAD,phiRAD,r,lmax,const,drtype,wbh,jflag)

% GRADPSHS determines the gradient of the potential pointwise at the
% location (r,\phi,\lambda). The calculation is done in the Earth-fixed
% frame with non-singular formulas (p.24 Mayer-Gürr,2006)
%
% IN:
%    field ... [n x m]   a priori gravity field in cs or sc format
%    lamRAD .. [t x 1]   longitude in [rad] 
%    phiRAD .. [t x 1]   latitude in [rad]
%    r ....... [t x 1]   radius in [m]. (r can be scalar)
%    lmax .... [1 x 1]   maximum degree 
%    const ... [2 x 1]   constants GM (element 1) and ae (element 2)
%    drtype .. [string]  defines the type of output
%                        - 'xyz'   are the derivatives towards {x,y,z}
%                                  in the cartesian frame
%                        - 'lscs'  are the derivatives towards {x,y,z}
%                                  but in the local spherical
%                                  coordinate frame (default) 
%    wbh ..... [1 x 1]   waitbar handle
%    jflag ... [1 x 1]   1 = reference field is subtracted (default)
%                        0 = reference field is not subtracted
% OUT:
%    gradV ... [t x 3]   gradient [V_x  V_y  V_z] in [m/s^2]
% 
% EXAMPLE:
%    load egm96; lmax = 200;
%    field =  vec2cs(lmax,EGM96(:,3),EGM96(:,4))
%    %% (extremly drifting) circular orbit:
%        x = (0:.01:10*pi)'; inc = -85;  omasc = x/15;
%        r = 7200000; z = 0*x;
%        X=multmatvek([cosd(z+inc), z, -sind(z+inc),z z+1,z, sind(z+inc), z, cosd(z+inc)],[r.*cos(x),r.*sin(x),z])
%        X=multmatvek([cosd(omasc), sind(omasc),z, -sind(omasc),  cosd(omasc), z, z, z, z+1],X)
%    [lamRAD,phiRAD,r] = cart2sph(X(:,1),X(:,2),X(:,3));
%    gradV_xyz = gradpshs(field,lamRAD,phiRAD,r,lmax,[],'xyz')
%    gradV_XYZ = gradpshs(field,lamRAD,phiRAD,r,lmax,[],'lscs')
%    figure
%    subplot(311);plot(gradV_XYZ); title('local frame'); legend('X','Y','Z')
%    subplot(312);plot(gradV_xyz); title('Earth-fixed frame'); legend('x','y','z')
%    subplot(313);plot(sum(gradV_XYZ.^2,2),'m'); hold on; plot(sum(gradV_xyz.^2,2),'k.');
%    title('absolut values in both frames'); legend('|gradV_{XYZ}|','|gradV_{xyz}|')
%
% USES: 
%    vec2cs, plm, normalklm, cs2sc, sc2cs, 
%    uberall/multmatvek, uberall/constants, uberall/checkcoor

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
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users)
%    2013-02-26: MR, doesn't calculate V_r V_phi V_lambda --> removed help
%                    text regarding this topic
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments
%    2013-01-23: MA, output in radian
%    2004-12-12: MW, change order of output arguments
%    2003-07-22: MW, initial version
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
if nargin < 5 || isempty(lmax),   lmax   = [];     end

% load necessary constants
if isempty(const)
    constants;
elseif numel(const) == 2
    GM = const(1);
    ae = const(2);
else
    error('CONST must be a 2x1 vector or empty.')
end


%% PREPARATION
% prepare the coordinates
[lamRAD,phiRAD,r] = checkcoor(lamRAD,phiRAD,r,ae,'pointwise');

% preparing the input field
[row, col] = size(field);
if (row~=col) && (col~=2*row-1)
   error('Input ''field'' not in cs or sc format');
elseif col==2*row-1
   field = sc2cs(field);
end
if isempty(lmax) || row-1 < lmax, lmax = row-1; end
field = field(1:lmax+1,1:lmax+1); 
field = cs2sc(field,0);
row   = size(field,1);

% substract reference field if requested
if jflag
    field = field - cs2sc(normalklm(lmax,'wgs84')); 
    warning('A reference field (WGS84) is removed from your coefficients')
end

% prepare the output
gradV = zeros(length(lamRAD),3);
didx  = 1:numel(lamRAD);
didx  = didx(all(~isnan([lamRAD phiRAD r]),2));

%% CALCULATION
maxlength = 512;
if ~isempty(wbh),
    wbhtext = get(get(get(wbh,'Children'),'Title'),'String');
    set(get(get(wbh,'Children'),'Title'),'String',[wbhtext ': ' num2str(0,'%02.2f') '%']);
    drawnow
end

% process data (piecewise)
parts = ceil(numel(didx)/maxlength);
for I = 1:parts
    % select the data
    % prepare cosine and sine --> cos(m*lam) and sin(m*lam) and  get
    % co-latitude 
    idx = didx((I-1)*maxlength+1:min(I*maxlength,numel(didx)));
    theRAD  = (pi/2 - phiRAD(idx));
    l   = (0:lmax);
    TF  = ae./r(idx) * ones(1,row);     % matrix of size(length(r),length(n))
    TF  = bsxfun(@power,TF,l+2);
    TSQ = realsqrt((2.*l'+1)./(2.*l'+3));
    TK  = GM./(2.*ae.^2);
    
    for m = 0:lmax
        l = (m:lmax)';
        if m == 0
            pClm  = field(:,row+m).*TSQ.*realsqrt((l+m+1).*(l+m+2).*2)  ;
            oClm  = field(:,row+m).*TSQ.*realsqrt((l-m+1).*(l+m+1))     ;
            pPQ   = TF.*(plm(l'+1,m+1,theRAD));
            oPQ   = TF.*(plm(l'+1,m,  theRAD));
            pcos  = cos((m+1).*lamRAD(idx));
            psin  = sin((m+1).*lamRAD(idx));
            ocos  = cos(m.*lamRAD(idx));
            % ----------------------------
            gradV(idx,1) = gradV(idx,1) -    TK.*(pPQ*pClm).*pcos;
            gradV(idx,2) = gradV(idx,2) -    TK.*(pPQ*pClm).*psin;
            gradV(idx,3) = gradV(idx,3) - 2.*TK.*(oPQ*oClm).*ocos;
            
        elseif m == 1
            pClm  = field(m+1:end,row+m).*TSQ(m+1:end).*realsqrt((l+m+1).*(l+m+2))     ;
            oClm  = field(m+1:end,row+m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l+m+1))     ;
            mClm  = field(m+1:end,row+m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l-m+2).*2)  ;
            pSlm  = field(m+1:end,row-m).*TSQ(m+1:end).*realsqrt((l+m+1).*(l+m+2))     ;
            oSlm  = field(m+1:end,row-m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l+m+1))     ;
            mSlm  = field(m+1:end,row-m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l-m+2).*2)  ;
            mPQ   = oPQ(:,2:end);
            oPQ   = pPQ(:,2:end);
            pPQ   = TF(:,2:end).*(plm(l'+1,m+1,theRAD));
            mcos  = ocos;
            ocos  = pcos;
            pcos  = cos((m+1).*lamRAD(idx));
            osin  = psin;
            psin  = sin((m+1).*lamRAD(idx));
            % ----------------------------
            gradV(idx,1) = gradV(idx,1) +    TK.*(mPQ*mClm).*mcos -    TK.*(pPQ*pClm).*pcos  -  TK.*(pPQ*pSlm).*psin;
            gradV(idx,2) = gradV(idx,2) -    TK.*(pPQ*pClm).*psin +    TK.*(mPQ*mSlm).*mcos  +  TK.*(pPQ*pSlm).*pcos;
            gradV(idx,3) = gradV(idx,3) - 2.*TK.*(oPQ*oClm).*ocos - 2.*TK.*(oPQ*oSlm).*osin;
            
        else
            pClm  = field(m+1:end,row+m).*TSQ(m+1:end).*realsqrt((l+m+1).*(l+m+2));
            oClm  = field(m+1:end,row+m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l+m+1));
            mClm  = field(m+1:end,row+m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l-m+2));
            pSlm  = field(m+1:end,row-m).*TSQ(m+1:end).*realsqrt((l+m+1).*(l+m+2));
            oSlm  = field(m+1:end,row-m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l+m+1));
            mSlm  = field(m+1:end,row-m).*TSQ(m+1:end).*realsqrt((l-m+1).*(l-m+2));
            mPQ   = oPQ(:,2:end);
            oPQ   = pPQ(:,2:end);
            pPQ   = TF(:,m+1:end).*(plm(l'+1,m+1,theRAD));
            mcos  = ocos;
            ocos  = pcos;
            pcos  = cos((m+1).*lamRAD(idx));
            msin  = osin;
            osin  = psin;
            psin  = sin((m+1).*lamRAD(idx));
            % ----------------------------
            gradV(idx,1) = gradV(idx,1) +    TK.*(mPQ*mClm).*mcos -    TK.*(pPQ*pClm).*pcos  +  TK.*(mPQ*mSlm).*msin  - TK.*(pPQ*pSlm).*psin;
            gradV(idx,2) = gradV(idx,2) -    TK.*(mPQ*mClm).*msin -    TK.*(pPQ*pClm).*psin  +  TK.*(mPQ*mSlm).*mcos  + TK.*(pPQ*pSlm).*pcos;
            gradV(idx,3) = gradV(idx,3) - 2.*TK.*(oPQ*oClm).*ocos - 2.*TK.*(oPQ*oSlm).*osin;
        end
    end
    
    % prepare output and rotate if desired
    if strcmpi(drtype,'lscs')
        cp = cos(phiRAD(idx));
        sp = sin(phiRAD(idx));
        cl = cos(lamRAD(idx));
        sl = sin(lamRAD(idx));
        gradV(idx,:) = multmatvek([cp.*cl,  cp.*sl,    sp,  -sp.*cl,  -sp.*sl,   cp,  -sl,  cl,  zeros(size(sp))],gradV(idx,:));
    end
    if any(size(drtype)==9)
        gradV(idx,:) = multmatvek(drtype(idx,:),gradV(idx,:));
    end

    % CleanUp
    clear idx l theRAD TF TSQ TK pClm oClm mClm pSlm oSlm mSlm mPQ oPQ pPQ mcos ocos pcos msin osin psin cp sp cl sl
    
    if ~isempty(wbh), 
        set(get(get(wbh,'Children'),'Title'),'String',[wbhtext ': ' num2str(I/parts*100,'%02.2f') '%']);
        drawnow
    end
    
end % end for I=1:parts
if ~isempty(wbh), set(get(get(wbh,'Children'),'Title'),'String',wbhtext); end

