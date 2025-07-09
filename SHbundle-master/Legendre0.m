function [PNM,dPNM] = Legendre0(lmax,normPnm)
% [PNM,dPNM] = Legendre0(lmax,normPnm)
%
% LEGENDRE0 calculates all Legendre functions and their first derivatives 
% at the equator in a triangle format. The zero elements of the odd functions
% are guaranteed by the modification of the recursion formula
%
% The normalization is by default the geodetic one; for the complex 
% inclination functions this can be changed by using 'normPnm'
%   'cplus':    Pnm_ = Pnm * (-1)^m / sqrt(2)  and Pn0_ = Pn0
%   'cminus':   Pnm_ = Pnm / sqrt(2)           and Pn0_ = Pn0
%   'cplusminus': combination of cplus and cminus in /S|C\ format
% This normalization is equivalent to the option 'geocomplex' in the
% routine LeNorm.
%
% IN:
%    lmax ........ maximum degree of the Legendre functions        [1 x 1]
%    normPnm...... optional change to complex normalization        [string]
%                   'cplus'|'cminus'|'cplusminus'
% 
% OUT:
%    P ........... normalized Legendre functions                   [N x N]                             
%    dP .......... first order derivatives of Legendre functions   [N x N]   
%
 
 
 
 
% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2016-05-23: MA, consider case lmax = 0
%    2016-05-18: MA, optional change of normalization
%    2016-05-10: MA, initial version
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
 
if lmax == 0
    PNM  = 1;
    dPNM = 0;
    return;
end
 
 
 
%% initialization for the "single point Legendre functions":
lmax1 = lmax + 1;
PNM(lmax1,lmax1)  = 0;
dPNM(lmax1,lmax1) = 0;
 
 
%% auxiliary values of the recursion
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
 
 
 
 
 
 
 
    %% == begin(LEGENDRE FUNCTIONS) =======================================   
 
    % sectorial functions:
    Pnnsint = HNM.*1.^(mvec-1);
    Pnn     = Pnnsint;  
    dPnn    = mvec.*0;  
 
    
    % definition of P{0,0} (derivative = 0)  
    PNM(1,1)   = Pnn(1);
       
    % calculation of P{1,0} and its derivative      
    PNM(2,1)   =  0;
    dPNM(2,1)  = -sqrt3;
    
    % calculation of P{1,1} and its derivative  
    PNM(2,2)   =  Pnn(2);
    dPNM(2,2)  =  dPnn(2);
 
    
    % initial values (degree-1) and (degree-2) for all orders m
    Pnm2        = PNM(1,:);
    Pnm1        = PNM(2,:);
    dPnm2       = dPNM(1,:);
    dPnm1       = dPNM(2,:);
 
    % coefficients (degree-1) 
    wnm1    = sqrt3;
    for n0  = 2:lmax
        % column and coefficients
        n1      = n0+1;
        wnm0    = WNM(n1,:);
        
        % recursion for tesseral and zonal derivatives
        dPnm0   = wnm0.*(-Pnm1-1./wnm1.*dPnm2);
        % sectorial derivative
        dPnm0(1,n1) = dPnn(n1);
        
        
        % recursion for tesseral and zonal functions
        Pnm0        = wnm0.*(-1./wnm1.*Pnm2);
 
        
        % correction of sectorial function
        Pnm0(1,n1)  = Pnn(n1);
 
        
        % update of the recursion
        Pnm2      = Pnm1;
        Pnm1      = Pnm0;
 
        dPnm2     = dPnm1;
        dPnm1     = dPnm0;
        wnm1      = wnm0;
       
        % matrix of Legendre functions and their derivatives
        PNM(n0+1,:) = Pnm1;
        dPNM(n0+1,:)= dPnm1;
    end
    %% == end(LEGENDRE FUNCTIONS) =========================================
    
if nargin > 1
    normPnm = lower(normPnm);
    if strcmp(normPnm(1),'c') == 1
        PNM(:,2:end) = PNM(:,2:end) * sqrt(1/2);
        dPNM(:,2:end)= dPNM(:,2:end) * sqrt(1/2);
    end
    switch normPnm
    case 'cplus'
        PNM(:,2:2:end) = (-1)*PNM(:,2:2:end);
        dPNM(:,2:2:end)= (-1)*dPNM(:,2:2:end);
    case 'cminus'
        PNM = fliplr(PNM);
        dPNM = fliplr(dPNM);
    case 'cplusminus'
        A = fliplr(PNM);
        B = fliplr(dPNM);
        PNM(:,2:2:lmax1) = (-1)*PNM(:,2:2:lmax1);
        dPNM(:,2:2:lmax1)= (-1)*dPNM(:,2:2:lmax1);
        
        PNM = [A,PNM(:,2:lmax1)];
        dPNM = [B,dPNM(:,2:lmax1)];
    otherwise
        warning('unknown normalization, please check options')
    end
        
 
end
 


