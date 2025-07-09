function [l,r] = degcorr(lmin,signal,varargin)

% DEGCORR determines the degree correlation between two or more models. The
% first input argument is the model to be tested against a number of other
% models. The output is the degree correlation for each combination. The
% length of each degree correlation depends on the input of length of the
% two models. 
%          
% IN:
%    lmin .... [1 x 1]   minimum degree to consider
%    signal .. [SC/CS]   set of SH coefficients, either in SC-triangle or 
%                        CS-square format
%    field1 .. [SC/CS]   analog  
%           
% OUT: 
%    l .... [cell]    degree vector for degree correlation
%    r .... [cell]    degree correlation
%
% EXAMPLE:
%    r = degcorr(lmin, signal, field1, field2,...)
%
% USES:
%    cs2sc, sc2cs

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), DoGE, UofC
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2004-09-28: MW, initial version
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

% ----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
% ----------------------------------------------------------------------------
% bring signal to SC format and set background values to zero if necessary
[row, col]=size(signal);
if (row~=col) && (col~=2*row-1)
    error('Input not in cs or sc format');
elseif col~=2*row-1
    signal = cs2sc(signal,0);
else
    signal = sc2cs(signal);
    signal = cs2sc(signal,0);
end
[row, ~] = size(signal);

% ----------------------------------------------------------------------------
% COMPUTATION
% ----------------------------------------------------------------------------
for k = 1:length(varargin)
    probant = signal;
    field   = varargin{k};
    
    % Perpare input data
    [rowf, colf] = size(field);
    if (rowf ~= colf) && (colf ~= 2 * rowf - 1)
        error('Input not in cs or sc format');
    elseif colf ~= 2 * rowf - 1
        field = cs2sc(field, 0);
    else
        field = sc2cs(field);
        field = cs2sc(field, 0);
    end
    [rowf, ~] = size(field);

    % cut to the same size
    if row <= rowf
        lmax    = row - 1;
        field   = field(1:lmax+1, rowf-lmax:2*lmax+rowf-lmax);
    elseif rowf < row
        lmax    = rowf - 1;
        probant = probant(1:lmax+1, row-lmax:2*lmax+row-lmax);
    end
    l{k} = lmin:lmax;

    % determine degree variance for the signal
    sigma_A = sqrt(sum(probant.^2, 2));
    sigma_A(sigma_A == 0) = 1e-30;

    % determine degree variance for the comparison
    sigma_B = sqrt(sum(field.^2, 2));
    sigma_B(sigma_B == 0) = 1e-30;

    % calculate degree correlation:
    erg  = sum(probant .* field,2) ./ (sigma_A .* sigma_B);
    r{k} = erg(lmin+1:end);    
end

