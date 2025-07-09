function varargout = clm2sc(clm, varargin)

% CLM2SC converts a list of coefficients in clm-format to /S|C\-format 
% and--if available--does the same with the standard deviations.
%
% IN:
%    clm ............ coefficients as matrix in the form [l m C S [sigma_C, sigma_S]]
%    'max_lm' ....... only coefficients up to a degree = order = max_lm are
%                     used. (default: read all coefficients (max_lm = inf)) 
%    'gcoef2sc' ..... provides functionality of deprecated function gcoef2sc, 
%                     see "OUT (gcoef2sc)" below. (default: false)
%
% OUT (standard):
%    sc ............. coefficients in /S|C\ format
%    maxGOout ....... maximum degree and order available
%    stdsc .......... standard deviations of the coefficients in /S|C\ (if
%                     available)
%
% OUT (gcoef2sc):
%    sc ............. coefficients in /S|C\ format
%    c .............. only C coefficients in lower triangular matrix
%    s .............. only S coefficients in lower triangular matrix
%
% EXAMPLE:
%    clm2sc(clm, 'max_lm', 30, 'gcoef2sc', false);

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias ROTH (MR), GIS, Uni Stuttgart 
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-10-22: MR, outsource parameter checking to uberall/GETOPT
%    2014-09-25: MR, merge behaviour of GCOEF2SC into ICGEM2SC and rename the 
%                    product into CLM2CS, update help text
%    2014-09-18: MR, add conversion for standard deviations of coefficients
%    2013-02-14: MR, modify the function that it works like promised in the
%                    help text, remove/replace now unnecessary comments
%    2012-06-26: MR, initial version
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

% define defaults and check optional parameters
defaultparams = {'max_lm', inf;
                 'gcoef2sc', false};
params = getopt(defaultparams, false, varargin);  

% scrap data of degree/order larger than params.max_lm
t = (clm(:, 1) > params.max_lm) | (clm(:, 2) > params.max_lm);
clm(t, :) = []; 

% cut params.max_lm down to the present maximum degree/order
t = max(max(clm(:, 1)), max(clm(:, 2)));
if t < params.max_lm
    params.max_lm = t;
end

% initialize output
if ~params.gcoef2sc
    varargout{1} = zeros(params.max_lm + 1, 2 * (params.max_lm + 1) - 1);
else % keep compatibility to the deprecated function params.gcoef2sc
    varargout{1} = ones(params.max_lm + 1, 2 * (params.max_lm + 1) - 1) * 1e-40;
    varargout{2} = ones(params.max_lm + 1, params.max_lm + 1) * 1e-40;
    varargout{3} = ones(params.max_lm + 1, params.max_lm + 1) * 1e-40;
end

% calculate indices
idx_s = sub2ind(size(varargout{1}), clm(:, 1) + 1, params.max_lm + 1 - clm(:, 2));
idx_c = sub2ind(size(varargout{1}), clm(:, 1) + 1, params.max_lm + 1 + clm(:, 2));

% put coeffs into matrix
varargout{1}(idx_s) = clm(:, 4); % ATTENTION! put first the S_lm, otherwise the S_l0 = 0 would overwrite C_l0!
varargout{1}(idx_c) = clm(:, 3);

if ~params.gcoef2sc
    varargout{2} = params.max_lm; 

    % if output argument stdsc is specified
    if (nargout == 3) && (size(clm, 2) == 6) % put standard deviations of coeffs into matrix
        % initialize output
        varargout{3} = zeros(params.max_lm + 1, 2 * (params.max_lm + 1) - 1);
        varargout{3}(idx_s) = clm(:, 6); % ATTENTION! put first the sigma_S_lm, otherwise the sigma_S_l0 = 0 would overwrite sigma_C_l0!
        varargout{3}(idx_c) = clm(:, 5);
    else
        varargout{3} = [];
    end
else % keep compatibility to the deprecated function params.gcoef2sc
    idx = sub2ind(size(varargout{2}), clm(:, 1) + 1, clm(:, 2) + 1);
    varargout{2}(idx)  = clm(:, 3);
    varargout{3}(idx)  = clm(:, 4);
    varargout{3}(:, 1) = 0; % S_l0 = 0 per definitionem
end
