function [outC,outS,nm] = cs2vec(in,ordflag)

% CS2VEC rearranges a field of spherical harmonic coefficients in 
% cs- or sc-format to a vector shape. The order will be:
%
% CASE I: degree-wise (default)
%        nm   = [0    1    1    2    2    2    3    3    3    3  ...
%                0    0    1    0    1    2    0    1    2    3  ...]
% CASE II: order-wise
%        nm   = [0    1    2    3    1    2    3    2    3    3  ...
%                0    0    0    0    1    1    1    2    2    3  ...]
%
% IN: 
%    in ....... [n,m]   input in cs- or sc-format
%    ordflag .. [bool]  false = degree-wise ordering (default)
%                       true  = order-wise ordering 
% 
% OUT: 
%    outC ..... [k,1]   coefficients in vector shape: cosine part
%    outS ..... [k,1]   coefficients in vector shape: sine part
%    nm ....... [2,k]   ordering vector
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
%    2012-03-12: MW, added ordering flag and completed help text
%    2007-05-02: MW, initial version
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

%% INPUT CHECK
% check if there is any NaN
narginchk(1,2);
if nargin < 2 || isempty(ordflag), ordflag = false;                       end
if any(isnan(in)), error('Input may not contain NaNs. Use Inf instead!'); end

% Size determination for field
[row, col]=size(in);
if (row~=col) && (col~=2*row-1)
   error('Input not in cs or sc format');
elseif col~=2*row-1
    % if data is in cs-format we transfer it to sc-format
   in = cs2sc(in,NaN);
   row = size(in,1);
else
   in = sc2cs(in);
   in = cs2sc(in,NaN);
end
lmax = row-1;

%% Processing
if ordflag
    % Order-wise ordering
    % -------------------------------------------
    % split sine and cosine part
    Cnm = in(:,lmax+1:end);
    Snm = [zeros(lmax+1,1) fliplr(in(:,1:lmax))];
    % -------------------------------------------
    % rearrange to a vector ordered by the degree
    outC = Cnm;
    outC = outC(:);
    outC(isnan(outC)) = [];
    outS = Snm;
    outS = outS(:);
    outS(isnan(outS)) = [];
    % -------------------------------------------
    % perpare the degree and order vector
    nelem = sum(0:lmax+1);
    nm    = zeros(2,nelem);
    sidx  = 1;
    for m = 0:lmax
        eidx            = sidx + lmax-m;
        nm(1,sidx:eidx) = m:lmax;
        nm(2,sidx:eidx) = m;
        sidx            = eidx + 1;
    end
    
else
    % Degree-wise ordering
    % -------------------------------------------
    % split sine and cosine part
    Cnm = in(:,lmax+1:end);
    Snm = [zeros(lmax+1,1) fliplr(in(:,1:lmax))];
    % -------------------------------------------
    % rearrange to a vector ordered by the degree
    outC = Cnm';
    outC = outC(:);
    outC(isnan(outC)) = [];
    outS = Snm';
    outS = outS(:);
    outS(isnan(outS)) = [];
    % -------------------------------------------
    % perpare the degree and order vector
    nc   = sum(1:lmax+1);  % number of elements
    i = 1;
    nm = zeros(2,nc);
    for l = 1:lmax+1
        for m = 1:l
            nm(1,i) = l-1;
            nm(2,i) = m-1;
            i = i+1;
        end
    end
end
outC = outC';
outS = outS';
