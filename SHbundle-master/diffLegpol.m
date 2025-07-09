function [pn, Pn,d1Pn,d2Pn,d3Pn] = diffLegpol(maxGrad,thetaRAD)

% DIFFLEGPOL calculates the non-normalized Legendre polynomials 
% P[n](t) = P[n,0](t) with t = cos(thetaRAD) including the derivatives w.r.t. 
% the argument t = cos(thetaRAD). 
%
% The routine uses the well-known recursion
%       (n+1)*P[n+1](t) = (2n+1)*t*P[n](t)-n*P[n-1](t)
% respectively
%       P[n+1](t) = ((2n+1)*t*P[n](t)-n*P[n-1](t))/(n+1).
%
% The polynomials are differentiated k = {1,2,3} times w.r.t. the argument 
% t and the 1. and 2. derivatives are tested by the Legendre differential 
% equation.
% 
% 
% IN:
%    maxGrad ... maximum degree of the Legendre polynomials                        [1x1]
%    thetaRAD .. co-latitude of the observation point                         [rad][1x1]
%
% OUT:
%    pn ........ Legendre polynomial of degree n = maxGrad
%    Pn ........ all calculated Legendre polynomials between degree 0 and degree n [NxN]
%    d1Pn ...... 1.derivatives of all Legendre polynomials                         [NxN]
%    d2Pn ...... 2. derivatives of all Legendre polynomials                        [NxN]
%                (checked with the Legendre-differential equation)
%    d3Pn ...... 3. derivatives of all Legendre polynomials                        [NxN]
%
% EXAMPLE:
%    clearvars  -except CONSTANT_EARTH_MODELL EXAMPLE_PATH
%    thetaRAD = 0:.01:pi; maxGrad = 10
%    pn0 = plm(maxGrad,0,thetaRAD);
%    [pn,Pn,dPn,ddPn] = diffLegpol(maxGrad,thetaRAD);
%    pn  = sqrt((2*maxGrad+1))*pn;
%    
%    figure;subplot(211);
%    plot(thetaRAD,pn0); hold on; plot(thetaRAD,pn,'r.');
%    legend('plm','diffLegpol')
%    subplot(212);
%    semilogy(thetaRAD,abs(pn-pn0),'k'); title('differences')
%
%    % Legendre-differential equation:
%    T = repmat(cos(thetaRAD),maxGrad+1,1);
%    Nvec = 0:maxGrad;
%    Nvec = Nvec.*(Nvec+1);
%    Nvec = diag(Nvec);
%    LDGL = Nvec * Pn-2*T.*dPn + (1-T.^2).*ddPn; 
%    figure
%    imagesc(thetaRAD*180/pi,0:maxGrad,LDGL); colorbar;
%    title('''Error'' of the Legendre differential equation')
%    xlabel('co-latitude [deg]'); ylabel('degree n')
%
% SEE ALSO:
%    LEGPOL, LEGENDREP, PLM
%

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%   2013-03-07: MA, last modification
%   2013-01-22: MA, initial version
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
% --------------------------------------------------------------------------
 
if nargin<2
    thetaRAD = (0:180)*pi/180;
end;
 
 
%% Check of Input
if max(abs(thetaRAD))>pi  || min(thetaRAD)<0
    warning('are the angels given in radiant?')
end;
if length(maxGrad)~=1 || rem(maxGrad,1)~=0 || maxGrad<0
    error('nonInteger:error','degree ''maxGrad'' must be integer, scalar and positive')
end;
 
thetaRAD = (thetaRAD(:))';
ltheta = length(thetaRAD);
t = cos(thetaRAD);
 
 
p0 = ones(1,ltheta);
p1 = t;
 
%% Initialization
Pn(maxGrad+1,ltheta)   = 0;
d1Pn(maxGrad+1,ltheta) = 0;
d2Pn(maxGrad+1,ltheta) = 0;
d3Pn(maxGrad+1,ltheta) = 0;
 
% Legendre polynomials of degree 0 and 1
Pn(1,:) = p0;
Pn(2,:) = p1;
 
%% Recursion formulas 
    
% Legendre polynomials of degree 2,3,4
Pn(3,:)  = 0.5*(3*t.^2-1);
p0  = 5/2*t.^3-3/2.*t;
p1  = 1/8*(35*t.^4-30*t.^2+3);
Pn(4,:) = p0;
Pn(5,:) = p1;

% 1. derivatives of Legendre polynomials of degree 1,2,3,4
dp0 = 15/2*t.^2-3/2;
dp1 = 35/2*t.^3 - 30/4*t; 
d1Pn(2,:) = 1;
d1Pn(3,:) = 3*t;
d1Pn(4,:) = dp0;
d1Pn(5,:) = dp1;

% 2. derivatives of Legendre polynomials of degree 2,3,4
ddp0 = 15*t;
ddp1 = 105/2*t.^2-30/4;

d2Pn(3,:) = 3;
d2Pn(4,:) = ddp0;
d2Pn(5,:) = ddp1;

% 3. derivatives of Legendre polynomials of degree 3,4
dddp0 = 15;
dddp1 = 105*t;
d3Pn(4,:) = dddp0;
d3Pn(5,:) = dddp1;

% higher orders:
for nn=4:maxGrad-1
    eta = (2*nn+1)/(nn+1);
    sigma = nn/(nn+1);
    p2 = eta.*t.*p1 -sigma*p0;

    dp2   = eta.*(t.*dp1 + p1) -sigma*dp0;
    ddp2  = eta.*(t.*ddp1 + 2*dp1) -sigma*ddp0;
    dddp2 = eta.*(t.*dddp1 + 3*ddp1) -sigma*dddp0;

    % Update 
    p0 = p1;
    p1 = p2;
    dp0 = dp1;
    dp1 = dp2;
    ddp0 = ddp1;
    ddp1 = ddp2;
    dddp0 = dddp1;
    dddp1 = dddp2;

    % Results
    Pn(nn+2,:) = p2;
    d1Pn(nn+2,:) = dp2;
    d2Pn(nn+2,:) = ddp2;
    d3Pn(nn+2,:) = dddp2;
end;
 
 
%% cut-off for a maximum degree n smaller than 4 
% (these values are created in the initialization already)    
Pn   = Pn(1:maxGrad+1,:);
d1Pn = d1Pn(1:maxGrad+1,:);
d2Pn = d2Pn(1:maxGrad+1,:);
d3Pn = d3Pn(1:maxGrad+1,:);
 
pn = Pn(maxGrad+1,:);
 


