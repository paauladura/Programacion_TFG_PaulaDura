function fhandle = sph_gradshs(GM, R, SHfield, lmax, fieldName,sat)
% fhandle = sph_gradshs(GM, R, SHfield, lmax, fieldName,sat)
%
% SPH_GRADSHS creates a SINGLE POINT HANDLE for calculating the gradient 
% of a potential described in spherical harmonics for a single point in
% space.
% 
% The potential V 
% 
%                &N     &n                                                                       & & 
% V = GM/R * SUM & SUM  &  (R/r)^{n+1} Pnm(cos(theta)) * [Cnm*cos(m*lambda) + Snm*sin(m*lambda)]}& & 
%                &n=0   &m=0                                                                     & &
% 
% is differentiated w.r.t spherical coordinates 
% 
%  /             d/dr            \
% |  -1/r            * d/dtheta   | {V}
%  \ 1/(r*sin(theta) * d/dlambda /
% 
% to achieve the gradient in a local north oriented frame, i.e.  
%       * the 1. axis is pointing in radial direction
%       * the 2. axis is oriented in towards the North pole/rotation axis and
%       * the 3. axis points tangential in lambda direction
%
% The gradient is transformed into the Earth-fixed frame by the rotation matrix
% 
%      | sintheta.*coslambda, -costheta.*coslambda, -sinlambda|
%      | sintheta.*sinlambda, -costheta.*sinlambda,  coslambda|
%      |            costheta,             sintheta,          0].
% 
% and can be rotated into the inertial frame by using GAST, xp,yp,...
% in the following.
%
% To determine the gradient, two steps are necessary: 
%      1. create the function handle (arbitrary name):
%           accEarth = sph_gradSHS(GM,R,SHfield, lmax,fieldName);
%      2. calculate the gradient
%           gradV = accEarth(lambdaRAD,thetaRAD,r)
%
%
% ==== ATTENTION ==========================================================
% !   The algorithm is optimized for the orbit integration and can
%     handle only single points in the orbit for the gradient. For further 
%     speed up, there is no test for r<R.
% !   The handle method is a newer feature of MATLAB, older versions or
%     alternatives like OCTAVE might be unable to run the code.
%
% IN (handle):
% GM:     mass(Earth) * gravitational constant                    [m^3/s^2][1x1]
% R:      radius of the spherical Earth                                 [m][1x1]
% SHfield: gravity field in spherical harmonic coefficients up to degree 
%          and order <lmax> in the quadratic [Cnm\Snm]-format           [-][NxN]
% lmax:   maximum degree of the Spherical Harmonic Synthesis (SHS)      [-][1x1]
% fieldName: name of the gravity field                                     ['string'] (optional)
%
% OUT 
% fhandle: function handle                                                 [@f]
% 
% -------------------------------------------------------------------------
%   IN: (gradient)
%   lambdaRAD:longitude of the observation point                      [rad][1x1]
%   thetaRAD: co-latitude of the observation point                    [rad][1x1]
%   r:        radius of the observation point                           [m][1x1]
%             (no check for r<R)
%
%   OUT:
%   gradV:  gradient of the field in the cartesian earth-fixed frame       [3x1]
% -------------------------------------------------------------------------
%
% sub-routines:
%    *   SHSgradient (intern)
%
% REMARKS:
%    *   The recursion is modified to determine [Pnm/sin(theta)] for m>0
%        to cirumvent singularities at the pole with a modified recusion 
%        (presented in HAINES)
%
% LITERATURE:
%        G.V. HAINES (1987)
%        "Computer Programs for Spherical Cap Harmonic Analysis of Potential
%        and General Fields"
%        Computer & Geosciences, Vol 14, No 4, pp. 413-447
%
% EXAMPLE:
%        clear all; clc;
%        R =  6378136; GM = 3.986004418e+14;
%        %% random field
%        lmax = 90;
%        field =  randn(lmax+1,lmax+1); field(1:2,1:2) = 0; field(1,1) = 1;
%        %% (extremly drifting) circular orbit:
%        x = (0:.01:pi)'; inc = -85;  omasc = x/15;
%        r = 7200000; z = 0*x; dataLength = numel(x), t = 1:dataLength;
%        X=multmatvek([cosd(z+inc), z, -sind(z+inc),z z+1,z, sind(z+inc), z, cosd(z+inc)],[r.*cos(x),r.*sin(x),z]);
%        X=multmatvek([cosd(omasc), sind(omasc),z, -sind(omasc), cosd(omasc), z, z, z, z+1],X);
%        [lambdaRAD,phiRAD,r] = cart2sph(X(:,1),X(:,2),X(:,3));
%        thetaRAD = pi/2-phiRAD;
%
%        % METHODE 1: pointwise implementation 
%        pwSHSgradient = sph_gradshs(GM,R,field, lmax,'random');
%        gradV_pkt = NaN(dataLength,3);
%        tic
%        for tt = 1:dataLength
%           gradV_pkt(tt,:) =  pwSHSgradient(lambdaRAD(tt),thetaRAD (tt),r(tt));
%        end  
%        sph_gradshsTOC = toc
%
%        % REFERENCE: gradpshs (pointwise used)
%        gradV = NaN(dataLength,3);
%        tic;  
%        for tt = 1:dataLength
%           gradV(tt,:) = gradpshs(field,lambdaRAD(tt),phiRAD(tt),r(tt),lmax,[GM,R],'xyz',[],0);    
%        end  
%        gradpshsTOC = toc;
%        subplot(211);
%        plot(t,gradV_pkt(:,1)); hold on;  plot(t,gradV,'.');
%        plot(t,gradV_pkt); legend('sp4SHSgradient.m','gradpshs.m')
%        title('Gradient in the Earth-fixed frame')
%        
%        subplot(212);
%        semilogy(t,abs(gradV-gradV_pkt));legend('"x"','"y"','"z"')
%        title('Difference between sph\_gradSHS.m and gradpshs.m')
%
% See also: GRADPSHS
%
% USES:
%    gradpshs, uberall\multmatvek

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2013-11-29: MA, add license information
%    2013-04-04: MA, speed-up, add help text
%    2012-02-12: MA, initial version
% ----------------------------------------------------------------------------
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

%% optional input arguments 
if nargin<6
    sat = 'single';
end
if nargin < 5
   fieldName = [];
end
 
%% check of dimensions
[row, col] = size(SHfield);
if row ~= col
    whos
    error('dimension:error', 'coefficient matrix ''SHfield'' must be in the quadratic [Cnm\\Snm]-format')
end
if nargin<4
    lmax = row-1;
end

% check gravity field
if row<9
    % increase quadatic field of SH-coefficients for n<8 (zero padding)
    % (to avoid dimension controll or errors in the Legendre functions)
    varhelp(9,9) = 0;
    varhelp(1:row,1:row) = SHfield;
    SHfield = varhelp;
    lmax = 8;
else
    lmax = min([row-1, lmax]);
    SHfield = SHfield(1:lmax+1, 1:lmax+1);
end

%% screen comment about the field and method
fprintf(1,'\ngradient is determined in the LNOF and transformed into Earth-fixed frame by rotation\n');
fprintf(1,'\nsingularities at the poles are overcome by a modified recursion for Pnm(cos(theta))/sin(theta)\n');
if ischar(fieldName)==0
    fprintf(1,'\ngravity field of the Earth up to degree and order %d\n', lmax);
else
    fprintf(1,'\ngravity field of the Earth: ''%s'' up to degree and order %d\n', fieldName,lmax);
end
lmax1 = lmax + 1;
 
%% initialization for the "single point Legendre functions":
RPNM(lmax1,lmax1)  = 0;
DRPNM(lmax1,lmax1) = 0;
mRPNM__sint(lmax1,lmax1) = 0;

%% auxilary values of the recursion
sqrt3 = sqrt(3);
mvec  = 0:lmax;

% coefficients of the 3-term recursion of Legendre
[m,n] = meshgrid(0:lmax);

DIAG    = diag(sqrt((2*(0:lmax)+1)./(2*(0:lmax))));
WNM     = tril(exp(0.5*(log(2*n+1)+log(2*n-1)-log(n+m)-log(n-m))),-1)+DIAG;
WNM     = WNM + triu(ones(lmax+1),1);
WNM(1,1)= 1;
WNM(2,1:2) = sqrt(3);
HNM     = cumprod(diag(WNM))';
clear DIAG n m

% scaling alread by GM/R
HNM     = GM/R*HNM;
GMsqrt3 = HNM(2);

%% separation of the SH-coefficients to 2 triangle matrices
if SHfield(1,1)~=1
    warning('unusual:parameter','C[0,0] is not equal to 1, possible reasons might a residual gravity field!')
end
Cnm     = tril(SHfield);
SHfield = triu(SHfield,1)';
Snm(:,2:lmax+1) = SHfield(:,1:lmax);
 
%% handle to calculate the gradient of the SH-expansion
%% selction between 1 or 2 satellites
switch lower(sat)
case 'single'
    fhandle = @SHSgradient1;
case 'double'
    RPNM = [RPNM,RPNM];
    DRPNM= [DRPNM,DRPNM];
    mRPNM__sint = [mRPNM__sint,mRPNM__sint];
    WNM = [WNM,WNM];
    ONE = ones(1,lmax1);
    fhandle = @SHSgradient2;
    ZweiMaxGrad = 2*lmax1;
    lmax2 =lmax+2;
otherwise
    error('nonImplemented:error','method is not implemented')
end
 
function [gradV] = SHSgradient1(lambdaRAD,thetaRAD,r,varargin)
% [gradV] = SHSgradient1(lambdaRAD,thetaRAD,r)
% 
% Single-point calculation of the gradient gradV of a potential field V.
%
% IN:
% lambdaRAD:longitude of the observation point                        [rad][1x1]
% thetaRAD: co-latitude of the observation point                      [rad][1x1]
% r:        radius of the observation point                             [m][1x1]
%
% OUT:
% gradV:  gradient of the field in the Cartesian Earth-fixed frame         [3x1]

    %% auxiliary values 
    costheta = cos(thetaRAD);
    sintheta = sin(thetaRAD);
    Rr       = R/r;
    RrCos    = Rr*costheta;
    RrSin    = Rr*sintheta;
    sinml    = sin(mvec'*lambdaRAD);    
    cosml    = cos(mvec'*lambdaRAD);
    Rr2      = Rr^2;
 
 
 
    %% == begin(LEGENDRE FUNCTIONS) =======================================   
    % recursion of the Legendre functions times the damping term 
    % Pnm(cos(theta)) * (R/r)^(n+1)
    
    % sectorial functions:
    Pnnsint = Rr.*HNM.*RrSin.^(mvec-1);
    Pnn     = Pnnsint.*RrSin;  
    dPnn    = mvec.*RrCos.*Pnnsint;  
    Pnnsint = Pnnsint*Rr;
    
    % definition of P{0,0} (derivative = 0)  
    RPNM(1,1)   = Pnn(1);
       
    % calculation of P{1,0} and its derivative      
    RPNM(2,1)   =  GMsqrt3*RrCos.*Rr;
    DRPNM(2,1)  = -GMsqrt3*RrSin.*Rr;
    
    % calculation of P{1,1} and its derivative  
    RPNM(2,2)   =  Pnn(2);
    DRPNM(2,2)  =  dPnn(2);
    mRPNM__sint(2,2) = Pnnsint(2);
    
    % initial values (degree-1) and (degree-2) for all orders m
    Pnm2        = RPNM(1,:);
    Pnm1        = RPNM(2,:);
    dPnm2       = DRPNM(1,:);
    dPnm1       = DRPNM(2,:);
    Pnm2sint    = mRPNM__sint(1,:);
    Pnm1sint    = mRPNM__sint(2,:);
    
    % coefficients (degree-1) 
    wnm1    = sqrt3;
    for n0  = 2:lmax
        % column and coefficients
        n1      = n0+1;
        wnm0    = WNM(n1,:);
        
        % recursion for tesseral and zonal derivatives
        dPnm0   = wnm0.*(RrCos.*dPnm1-RrSin.*Pnm1-1./wnm1.*dPnm2.*Rr2);
        % sectorial derivative
        dPnm0(1,n1) = dPnn(n1);
        
        
        % recursion for tesseral and zonal functions
        Pnm0        = wnm0.*(RrCos.*Pnm1-1./wnm1.*Pnm2.*Rr2);
        Pnm0sint    = wnm0.*(RrCos.*Pnm1sint-1./wnm1.*Pnm2sint.*Rr2);
        
        % correction of sectorial function
        Pnm0(1,n1)  = Pnn(n1);
        Pnm0sint(1,n1) = Pnnsint(n1);
        
        % update of the recursion
        Pnm2      = Pnm1;
        Pnm1      = Pnm0;
        Pnm2sint  = Pnm1sint;
        Pnm1sint  = Pnm0sint;
        dPnm2     = dPnm1;
        dPnm1     = dPnm0;
        wnm1      = wnm0;
       
        % matrix of Legendre functions and their derivatives
        RPNM(n0+1,:) = Pnm1;
        mRPNM__sint(n0+1,:) = mvec.*Pnm1sint;
        DRPNM(n0+1,:)= dPnm1;
    end
    %% == end(LEGENDRE FUNCTIONS) =========================================
 
 
    %% gradient in the LNOF
    gradV(3,1) =  sum((-mRPNM__sint.*Cnm)*sinml + (mRPNM__sint.*Snm)*cosml);
    gradV(2,1) = -sum((DRPNM.*Cnm)*cosml + (DRPNM.*Snm)*sinml) ;
    gradV(1,1) = -(mvec+1)*((RPNM.*Cnm)*cosml + (RPNM.*Snm)*sinml) ;
 
    %% transformation to Cartesian coordinates
    coslambda = cosml(2);
    sinlambda = sinml(2);
    Rl2e = [ sintheta.*coslambda, -costheta.*coslambda, -sinlambda;
            sintheta.*sinlambda, -costheta.*sinlambda, coslambda;
            costheta, sintheta, 0];
 
   
    gradV = Rl2e*gradV;   
 
    % scaling 
    gradV = gradV/r;   
    
 
end


function [gradV] = SHSgradient2(lambdaRAD,thetaRAD,r,varargin)
% [gradV] = SHSgradient2(lambdaRAD,thetaRAD,r,varargin)
% 
% Double-point calculation of the gradient gradV of a potential field V.
%
%
% IN:
% lambdaRAD:longitude of the observation point                        [rad][2x1]
% thetaRAD: co-latitude of the observation point                      [rad][2x1]
% r:        radius of the observation point                             [m][2x1]
%
% OUT:
% gradV:  gradient of the field in the cartesian earth-fixed frame         [3x2]






    %% auxiliary values 
    costheta = cos(thetaRAD);
    sintheta = sin(thetaRAD);
    coslambda= cos(lambdaRAD);
    sinlambda= sin(lambdaRAD); 
    Rr       = R./r;
    RrCos    = Rr.*costheta;
    RrSin    = Rr.*sintheta;
    sinml(:,1)    = sin(mvec'*lambdaRAD(1));    
    cosml(:,1)    = cos(mvec'*lambdaRAD(1));
    sinml(:,2)    = sin(mvec'*lambdaRAD(2));    
    cosml(:,2)    = cos(mvec'*lambdaRAD(2));
 


    %% == begin(LEGENDRE FUNCTIONS) =======================================  
    % recursion of the Legendre functions times the damping term 
    % Pnm(cos(theta)) * (R/r)^(n+1)
    
    % sectorial functions:
    Pnnsint(1,:) = Rr(1).*HNM.*RrSin(1).^(mvec-1);
    Pnn(1,:)   = Pnnsint(1,:).*RrSin(1);  
    dPnn(1,:)   = mvec.*RrCos(1).*Pnnsint(1,:);  
    Pnnsint(1,:) = Pnnsint(1,:)*Rr(1);
    
    Pnnsint(2,:) = Rr(2).*HNM.*RrSin(2).^(mvec-1);
    Pnn(2,:)   = Pnnsint(2,:).*RrSin(2);  
    dPnn(2,:)   = mvec.*RrCos(2).*Pnnsint(2,:);  
    Pnnsint(2,:) = Pnnsint(2,:)*Rr(2);
    
    
    
    % definition of P{0,0} (derivative = 0)  
    RPNM(1,1)    = Pnn(1,1);
    RPNM(1,lmax1+1)    = Pnn(2,1);   
    % calculation of P{1,0} and its derivative   
    RPNM(2,1) =  GMsqrt3*RrCos(1).*Rr(1);
    DRPNM(2,1)= -GMsqrt3*RrSin(1).*Rr(1);
    
    RPNM(2,lmax2) =  GMsqrt3*RrCos(2).*Rr(2);
    DRPNM(2,lmax2)= -GMsqrt3*RrSin(2).*Rr(2);
    
    % calculation of P{1,1} and its derivative  
    RPNM(2,2) =  Pnn(1,2);
    DRPNM(2,2) = dPnn(1,2);
    mRPNM__sint(2,2) = Pnnsint(1,2);

    RPNM(2,lmax1+2) =  Pnn(2,2);
    DRPNM(2,lmax1+2) = dPnn(2,2);
    mRPNM__sint(2,lmax1+2) = Pnnsint(2,2);
    
    % initial values (degree-1) and (degree-2) for all orders m
    Pnm2        = RPNM(1,:);
    Pnm1        = RPNM(2,:);
    dPnm2       = DRPNM(1,:);
    dPnm1       = DRPNM(2,:);
    Pnm2sint    = mRPNM__sint(1,:);
    Pnm1sint    = mRPNM__sint(2,:);
    
    % coefficients (degree-1) 
    wnm1    = sqrt3;
    RrSin = [ONE*RrSin(1),ONE*RrSin(2)];
    RrCos = [ONE*RrCos(1),ONE*RrCos(2)];
    Rr = [ONE*Rr(1),ONE*Rr(2)];
    Rr2 = Rr.^2;
    for n0  = 2:lmax
        % column and coefficients
        n1      = n0+1;
        wnm0    = WNM(n1,:);
        
        % recursion for tesseral and zonal derivatives

        dPnm0   = wnm0.*(RrCos.*dPnm1-RrSin.*Pnm1-1./wnm1.*dPnm2.*Rr2);
        % sectorial derivative
        dPnm0(1,n1) = dPnn(1,n1);
        dPnm0(1,lmax1+n1) = dPnn(2,n1);
        
        % recursion for tesseral and zonal functions
        Pnm0    = wnm0.*(RrCos.*Pnm1-1./wnm1.*Pnm2.*Rr2);
        Pnm0sint    = wnm0.*(RrCos.*Pnm1sint-1./wnm1.*Pnm2sint.*Rr2);
        
        % correction of sectorial function
        Pnm0(1,n1) = Pnn(1,n1);
        Pnm0sint(1,n1) = Pnnsint(1,n1);
        
        Pnm0(1,lmax1+n1) = Pnn(2,n1);
        Pnm0sint(1,lmax1+n1) = Pnnsint(2,n1);
        
        
        % update of the recursion
        Pnm2    = Pnm1;
        Pnm1    = Pnm0;
        Pnm2sint= Pnm1sint;
        Pnm1sint= Pnm0sint;
        dPnm2   = dPnm1;
        dPnm1   = dPnm0;
        wnm1    = wnm0;
       
        % matrix of Legendre functions and their derivatives
        RPNM(n0+1,:) = Pnm1;
        mRPNM__sint(n0+1,:) = [mvec,mvec].*Pnm1sint;
        DRPNM(n0+1,:)= dPnm1;
    end
    %% == end(LEGENDRE FUNCTIONS) =========================================

 


    %% gradient in the LNOF
    gradV(3,1) =  sum((-mRPNM__sint(:,1:lmax1).*Cnm)*sinml(:,1) + (mRPNM__sint(:,1:lmax1).*Snm)*cosml(:,1));
    gradV(2,1) = -sum((DRPNM(:,1:lmax1).*Cnm)*cosml(:,1) + (DRPNM(:,1:lmax1).*Snm)*sinml(:,1)) ;
    gradV(1,1) = -(mvec+1)*((RPNM(:,1:lmax1).*Cnm)*cosml(:,1) + (RPNM(:,1:lmax1).*Snm)*sinml(:,1)) ;

    gradV(3,2) =  sum((-mRPNM__sint(:,lmax2:ZweiMaxGrad).*Cnm)*sinml(:,2) + (mRPNM__sint(:,lmax2:ZweiMaxGrad).*Snm)*cosml(:,2));
    gradV(2,2) = -sum((DRPNM(:,lmax2:ZweiMaxGrad).*Cnm)*cosml(:,2) + (DRPNM(:,lmax2:ZweiMaxGrad).*Snm)*sinml(:,2)) ;
    gradV(1,2) = -(mvec+1)*((RPNM(:,lmax2:ZweiMaxGrad).*Cnm)*cosml(:,2) + (RPNM(:,lmax2:ZweiMaxGrad).*Snm)*sinml(:,2)) ;
    
    
    %% transformation to cartesian coordinates
    Rl2e = [ sintheta(1).*coslambda(1), -costheta(1).*coslambda(1), -sinlambda(1);
            sintheta(1).*sinlambda(1), -costheta(1).*sinlambda(1), coslambda(1);
            costheta(1), sintheta(1), 0];

   
    gradV(:,1) = Rl2e*gradV(:,1)/r(1);   
    
        Rl2e = [ sintheta(2).*coslambda(2), -costheta(2).*coslambda(2), -sinlambda(2);
            sintheta(2).*sinlambda(2), -costheta(2).*sinlambda(2), coslambda(2);
            costheta(2), sintheta(2), 0];
    
    
    gradV(:,2) = Rl2e*gradV(:,2)/r(2);   
  
    

  end


end
