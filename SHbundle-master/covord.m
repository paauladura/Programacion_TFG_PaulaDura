function covcell = covord(fname, pth, lmax, vartype)

% COVORD converts a covariance matrix of a spherical harmonic degree
% expansion into a structure of four variables with each variable 
% containing cells. The number cells equal the numbers of the orders.
% The input covariance matrix must be in order-leading ordering format.
%
% IN:
%    fname .... File name of the corresponding covariance matrix
%    pth ...... Path of the file location
%    lmax ..... Maximum degree of spherical harmonic expansion.
%    vartype .. Variance type to be extracted: 0=full; 1=block-diagonal;
%               2=diagonal;
%
% OUT: 
%    covcell .. Structure containing three variables (cc, cs, sc, ss), 
%               with containing cells 
%

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2008-01-15:  BD, initial version
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


if nargin > 4,  error('Too many input arguments'),  end
if nargin < 3,  error('Maximum spherical harmonic degree of expansion missing'), end
if nargin==3,   vartype = 0; end

if ~ischar(fname) 
    error('Filename is not a character array')
end

if ~isempty(pth) && ~ischar(pth), 
    error('Path name is neither empty nor a string.')
end
if mod((lmax*10),10) ~= 0
    display('Warning: Maximum degree of the spherical harmonic expansion is not an integer. \n It will be rounded to the nearest integer')
    lmax = floor(lmax);
end

vcm = load([pth,fname]);
fld = fieldnames(vcm);

if size(vcm.(fld{1,1}),1) ~= (lmax+1)^2
    error('Matrix size does not match the maximum degree of spherical harmonic expansion')
end

% vcm.(fld{1,1}) = tril(vcm.(fld{1,1}));

if vartype == 0 % Full matrix extracted
    covcell = struct('cc',0,'sc',0,'ss',0,'cs',0);
    covcell.cc = cell(lmax+1);
    covcell.cs = cell(lmax+1);
    covcell.sc = cell(lmax+1);
    covcell.ss = cell(lmax+1);

    cindx2 = cumsum(fliplr(1:lmax+1))';
    cindx1 = [1;cindx2(1:end-1)+1];
    sindx2 = cumsum(fliplr(1:lmax))' + cindx2(end);
    sindx1 = [cindx2(end)+1;sindx2(1:end-1)+1];

    covcell.cc(:,1) = mat2cell(full(vcm.(fld{1,1})(cindx1(1):cindx2(end),cindx1(1):cindx2(1))), fliplr(1:lmax+1)',lmax+1);
    
    covcell.cs(:,1) = mat2cell(zeros(sum(1:lmax+1),lmax+1),(lmax+1:-1:1),lmax+1);
    covcell.cs(1,2:end) = mat2cell(full(vcm.(fld{1,1})(cindx1(1):cindx2(1),sindx1(1):sindx2(end))), lmax+1,fliplr(1:lmax)');
    
    covcell.sc(1,:) = mat2cell(zeros(lmax+1,sum(1:lmax+1)),lmax+1,(lmax+1:-1:1));
    covcell.sc(2:end,1) = mat2cell(full(vcm.(fld{1,1})(sindx1(1):sindx2(end),cindx1(1):cindx2(1))), fliplr(1:lmax)',lmax+1);

    covcell.ss(:,1) = mat2cell(zeros(sum(1:lmax+1),lmax+1),(lmax+1:-1:1)',lmax+1);
    covcell.ss(1,:) = mat2cell(zeros(lmax+1,sum(1:lmax+1)),lmax+1,(lmax+1:-1:1));
    
    for i = 1:lmax
        covcell.cc(1:end,i+1) = mat2cell(full(vcm.(fld{1,1})(cindx1(1):cindx2(end),cindx1(i+1):cindx2(i+1))), fliplr(1:lmax+1)',lmax+1-i);
        covcell.cs(2:end,i+1) = mat2cell(full(vcm.(fld{1,1})(cindx1(2):cindx2(end),sindx1(i):sindx2(i))), fliplr(1:lmax)',lmax+1-i);
        covcell.sc(2:end,i+1) = mat2cell(full(vcm.(fld{1,1})(sindx1(1):sindx2(end),cindx1(i+1):cindx2(i+1))), fliplr(1:lmax)',lmax+1-i);
        covcell.ss(2:end,i+1) = mat2cell(full(vcm.(fld{1,1})(sindx1(1):sindx2(end),sindx1(i):sindx2(i))), fliplr(1:lmax),lmax+1-i);
    end
elseif vartype == 1  % Only block-diagonal extracted
    covcell = struct('cc',0,'ss',0);
    covcell.cc = cell(lmax+1,1);
    covcell.ss = cell(lmax+1,1);

    cindx2 = cumsum(fliplr(1:lmax+1))';
    cindx1 = [1;cindx2(1:end-1)+1];
    sindx2 = cumsum(fliplr(1:lmax))' + cindx2(end);
    sindx1 = [cindx2(end)+1;sindx2(1:end-1)+1];

    covcell.cc{1} = vcm.(fld{1,1})(cindx1(1):cindx2(1),cindx1(1):cindx2(1));
    covcell.ss{1} = sparse(zeros(lmax+1));
    for i = 1:lmax
        covcell.cc{i+1} = vcm.(fld{1,1})(cindx1(i+1):cindx2(i+1),cindx1(i+1):cindx2(i+1));
        covcell.ss{i+1} = vcm.(fld{1,1})(sindx1(i):sindx2(i),sindx1(i):sindx2(i));
    end
elseif vartype == 2 % Only diagonal elements extracted
    covcell = struct('cc',0,'ss',0);
    covcell.cc = cell(lmax+1,1);
    covcell.ss = cell(lmax+1,1);
    
    indx = sum(1:lmax+1);
    
    vcm = full(diag(vcm.(fld{1,1})));
    covcell.cc(:,1) = mat2cell(vcm(1:indx),fliplr(1:lmax+1)',1);
    covcell.ss(2:end,1) = mat2cell(vcm(indx+1:end),fliplr(1:lmax)',1);
    covcell.ss{1} = zeros(lmax+1,1);
end
