function [longE, idx] = long360toeast(long360, varargin)

% LONG360TOEAST converts longitude values in 0 to 360 to -180 to 180.
%
% longE = long360toeast(long360);
% longE = long360toeast(long360, OPTIONS,VALUE);
%
% INPUT
% long360   -   Longitude values in radians [n x 1] or [n x m]
% 
% OPTIONS       VALUES
% -------       ------
% 'sort'        true false[def.]
% 'dim'         Dimension to sort -- 1 [def.] [integer]
% 'mode'        'ascend'[def.] 'descend'
%
%
% OUTPUT
% longE     -   Longitude values in radians
% idx       -   Sorting index
%
%------------------------------------------------------------------------------
%
%

% Balaji Devaraju. 8 May 2015. IfE Hannover.

%------------------------------------------------------------------------------

if nargin < 1
    error('Insufficient input arguments')
end

def = { 'sort', false;
        'dim', 1;
        'mode', 'ascend'};

params = getopt(def, false, varargin);

tmp     = long360;
idx     = (tmp > pi);
if ~isempty(idx)
    tmp(idx) = tmp(idx) - 2*pi;
end
idx     = [];

if params.sort
    [longE, idx] = sort(tmp, params.dim, params.mode);
else
    longE = tmp;
end
