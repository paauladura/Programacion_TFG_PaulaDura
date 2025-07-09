function T = updwcont(T,r,potTr,potTrr,potTrrr,potTrrrr,r_mean)

% function to do upward and downward continuation to a sphere with mean radius 
%
% IN:
%    T ......... disturbing orbit along the orbit in [m^2/s^2]                   [m x 1]   
%    r ......... radius [m]                                                      [m x 1]   
%    potTr ..... first derivative of the potential from a model in [m/s^2]       [m x 1]   
%    potTrr .... second derivative of the potential from a model in [1/s^2]      [m x 1]   
%    potTrrr ... third derivative of the potential from a model in [1/ms^2]      [m x 1]     
%    r_mean .... mean radius [m]                                                 [1 x 1]   
% 
% OUT:
%    T ......... disturbing orbit at mean height in [m^2/s^2]                    [m x 1]   

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), DoGE, UofC
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2004-05-08: MW, initial version
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
% INPUT CHECK
% -------------------------------------------------------------------------
if nargin < 7, r_mean = potTrrrr; end

% -------------------------------------------------------------------------
% CONTINUATION
% -------------------------------------------------------------------------
dr     = r_mean-r;
T      = T  +  potTr.*dr  +  1/2.*potTrr.*dr.^2   +   1/6.*potTrrr.*dr.^3;
if nargin == 7, T = T + 1/24.*potTrrrr.*dr.^4; end

% That's it
