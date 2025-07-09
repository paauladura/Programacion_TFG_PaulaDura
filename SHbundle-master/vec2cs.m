function out = vec2cs(nm,inC,inS)

% VEC2CS rearranges a vector shaped set of spherical harmonic
% coefficients into cs-format. The vector order is given by the ordering
% vector.
%
% IN:
%    nm ...... ordering vector                                      [2,k]  
%              maximum degree - degree-wise ordering is assumed     [1,1]  
%    inC ..... coefficients in vector shape: cosine part            [1,k]  
%    inS ..... coefficients in vector shape: sine part              [1,k]  
% 
% OUTPUT: 
%    out ..... output in cs- or sc- format                          [n,m]  
%
% USES:
%    sc2cs
%
% SEE ALSO:
%    cs2vec

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), DoGE, UofC
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-03-12: MW, change input lmax to nm-vector
%    2007-05-02: MW, initial version
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

%% Input check
narginchk(3,3);
if size(nm,1)~=2 && ~isscalar(nm), error('nm must be a 2 x k vector or a scalar.');    end 
if ~isscalar(nm),
    if numel(inC)~=size(nm,2), error('Dimension of nm and inC do not match.');         end
    if numel(inS)~=size(nm,2), error('Dimension of nm and inS do not match.');         end
else
    lmax = nm;
    nm   = sum(1:lmax+1);  % number of elements
    if numel(inC)~=nm, error('Dimension of nm and inC do not match.');         end
    if numel(inS)~=nm, error('Dimension of nm and inS do not match.');         end
end

%% Preparation
% perpare the degree and order vector
if isscalar(nm),
    idx  = 1;
    nm   = zeros(2,nm);
    for l = 1:lmax+1
        for m = 1:l
            nm(1,idx) = l-1;
            nm(2,idx) = m-1;
            idx = idx+1;
        end
    end
end

%% Sort coefficients
lmax = max(nm(1,:));
outC = zeros(lmax+1);
outS = zeros(lmax+1);
for idx = 1:size(nm,2);
    outC(nm(1,idx)+1,nm(2,idx)+1) = inC(idx);
    outS(nm(1,idx)+1,nm(2,idx)+1) = inS(idx);
end
out = sc2cs([fliplr(outS(:,2:end)) outC]);

