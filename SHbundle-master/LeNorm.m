function [hnm,success]=LeNorm(order,normPnm)
 
% LENORM determines the order-dependent factor between the 
% Legendre functions in the GEODETIC normalization
% 
%       Pnm_geo = norm_geo * pnm(t)                   [1]
% 
% and some alternavtive conventions
% 
%       Pnm_alt = norm_alt * pnm(t)                   [2]
% 
% where pnm(t) are  the NON-NORMALIZED Legendre functions:
%
% pnm(t):= 1/(2^n*n!)*(1-t^2)^(m/2) * diff((t^2-1)^n, t, n+m)
% 
% Only the factor between [1] and [2] is determined. Hence, the conversion is
% 
%        Pnm_alt = hnm * Pnm_geo
%
% ====================================================================
% The geodetic normalization 
%
%                / sqrt(2 * (2n+1) * (n-m)!/(n+m)!)   for m~=0
%   Nnm_geo    ={ 
%                \ sqrt(2n+1)                         if m==0
% 
% leads to the integral relation (of equal order and degree)
%
%     1/                              / 4 if m ~= 0
%     | Pnm_geo(t) * Pnm_geo(t) dt = {
%     /-1                             \ 2 if m  = 0
%
% The functions of the SHBUNDLE and all spherical harmonic coefficients of
% the gravity field models require the geodetic normalization. 
% 
% In the complex representation some features are easier to implement, e.g. 
% the rotation by the SO(3) coefficients (= Wigner-functions * exponentials). 
%
% ====================================================================
% There are several alternative conventions 
%
% A) MATHEMATICA (for complex representation):
% 
%              / (-1)^m*sqrt((2n+1)/(4*pi) * (n-m)!/(n+m)!)  for m>0
%   Nnm_com = {
%              \        sqrt((2n+1)/(4*pi) * (n-m)!/(n+m)!)  for m<0.
% 
%   => hnm = { (-1)^m * 1/(2*sqrt(2*pi)) | 1/2*sqrt(pi) | 1/(2*sqrt(2*pi)) } 
%                     m > 0              | m = 0        |       m< 0   
%
% B) geocomplex (complex standard in geodetic literature, cf. RUMMEL's lecture notes)
% 
%              / (-1)^m*sqrt((2n+1) * (n-m)!/(n+m)!)  for m>=0
%   Nnm_com = {
%              \        sqrt((2n+1) * (n-m)!/(n+m)!)  for m<0.
%
%   => hnm =      {  (-1)^m * 1/sqrt(2)  |  1           |  1/sqrt(2)  }
%                     m > 0              | m = 0        |       m< 0   
% C) ylmsimons (real representation; articles of Slepian basis functions)
%
%                                  / sqrt(2 * (2n+1) * (n-m)!/(n+m)!)   for m~=0
%   Nnm_com    =(-1)^m/sqrt(4*pi) { 
%                                  \ sqrt(2n+1)                         if m==0
%
%   => hnm =      (-1)^m/sqrt(4*pi)
% (F. Simons splits the normalization, the factor sqrt(2) is in Ynm, but not in Pnm)
%   
% IN:
%    order ........ order of the associated Legendre functions Pnm of degree n  [NxM]
%                   ( -n <= order<= n ) 
%    normPnm ....... select the alternative normalization                       ['string']
%                    * 'mathematica' (default due to prior project)
%                    * 'geocomplex'
%                    * 'ylmsimons'

% OUT:
%    hnm .......... factor for: Pnm_alt = Pnm_geo * LeNorm(order)               [Nx1]
%    success ...... returns 'false', if the normalization is not implemented
%
% REMARKS:
%    Please note, that there are several versions of complex normalization
%    which differ in the factors (-1)^m and/or 1/(4*pi). The default version 
%    here is identical with the one in MATHEMATICA.
% 
%
% EXAMPLE:
%    thetaRAD = 0:.01:pi; maxGrad = 10; normPnm = 'mathematica'
%    [pnm_real,delVar,delVar,DegreeOrder] = legendreP(maxGrad,thetaRAD);
%    pnm_comp_plus = bsxfun(@times,pnm_real,LeNorm(DegreeOrder(:,2),normPnm))
%    pnm_comp_minus = bsxfun(@times,pnm_real,LeNorm(-DegreeOrder(:,2),normPnm))

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2021-10-20: MA, new output 'success', default values hnm = 1
%    2015-07-07: MA, included 2 complex normalization, brushed up help text
%    2013-01-22: MA, translation of help text
%    2010-08-12: MA, factor 4*pi
%    2007-01-30: MA, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
 
if nargin < 2
    normPnm = 'mathematica';
end
 
if any(rem(order,1) ~= 0) == 1
    error('nonInteger:error','order must be integer')
end
 
success = true;
switch lower(normPnm)
case 'mathematica'
    % normalization in MATHEMATICA
    sqrtpi = 1/(2*sqrt(2*pi));
    
    % negative orders:
    hnm  = ones(size(order))*sqrtpi;
    
    % order zero
    ii = order==0;
    hnm(ii) = sqrt(2)*hnm(ii);
    
    % positive orders:
    ii = order>0;
    hnm(ii) = (-1).^order(ii)* sqrtpi;
case 'geocomplex'
    % complex normalization in geodesy 
    % negative orders:
    hnm  = ones(size(order));
    
    ii = order~=0;
    hnm(ii) = hnm(ii)/sqrt(2);
    
    % positive orders:
    ii = order>0;
    hnm(ii) = (-1).^order(ii).*hnm(ii);
case 'ylmsimons'
    % normalization of spherical harmonic functions in SIMON's articles about Slepian basis functions
    hnm  = (-1).^order / sqrt(4*pi); 
otherwise
    success = false;
    fprintf(1,'LENORM.M: unknown option <<%s>>: normalization remains unchanged',normPnm)
    hnm = ones(size(order));
end
 
 
 


