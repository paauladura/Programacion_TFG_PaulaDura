function [ylmc,ylms] = ylm(l,m,thetaRAD,lambdaRAD)

% YLM returns real normalized spherical harmonic.
%
% USE [ylmc,ylms] = ylm(l,m)                calculates on a 5x5 degree grid
%     [ylmc,ylms] = ylm(l,m,thetaRAD,lambdaRAD)  calculates on the prescribed grid
%
% IN:
%    l ......... spherical harmonic degree  (scalar)         
%    m ......... order                      (scalar)
%    thetaRAD .. co-latitude [radian]       (vector)
%    lambdaRAD . longitude   [radian]       (vector)
%
% OUT:
%    ylmc ...... cosine part
%    ylms ...... sine part
% 
% USES:
%    plm
%
% SEE ALSO:
%    PLM

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAGP, TU Munich
%    Matthias WEIGELT (MW),  GI, Uni Stuttgart 
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-29: MA, comments
%    2013-01-23: MA, input in radian
%    1997-02-21: NS, brushing up and redefinition of checks
%    1993-??-??: NS, initial version
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

% Check out L and M
if max(size(l)) > 1     , error('L must be scalar')              , end
if max(size(m)) > 1     , error('M must be scalar')              , end
if l < 0                , error('L < 0 ???')                     , end
if m < 0                , error('M < 0 ???')                     , end
if m > l                , error('M > L ???')                     , end

% Check out THETARAD and LAMBDARAD
if nargin == 2
   thetaRAD  = linspace(0,pi,37);				% create 5x5 grid
   lambdaRAD = linspace(0,2*pi,73);
end
if min(size(thetaRAD)) > 1 , error('''thetaRAD'' must be vector')              , end
if min(size(lambdaRAD)) > 1, error('''lambdaRAD'' must be vector')               , end
if min(thetaRAD)  <   0    , error('''thetaRAD'' < 0')                         , end
if max(thetaRAD)  > pi    ,  error('''thetaRAD'' must be given in radian')     , end
if max(lambdaRAD)  > 2*pi    ,  error('''lambdaRAD'' must be given in radian')   , end
% Convert THETARAD into standing vector and LAMBDARAD into lying one.
thetaRAD  = thetaRAD(:);				% column vector now
lambdaRAD = lambdaRAD(:)';			    % row vector now

% Actual computations
cosml = cos(m*lambdaRAD);
sinml = sin(m*lambdaRAD);
p     = plm(l,m,thetaRAD);

ylmc = p * cosml;
ylms = p * sinml;
