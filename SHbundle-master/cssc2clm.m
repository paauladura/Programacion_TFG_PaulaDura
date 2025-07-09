function vec = cssc2clm(mat,lmax)

% CSSC2CLM converts CS, SC format matrices and degree variance information
% to Colombo ordering vectors [l m Clm Slm]
%
% IN:
%    mat ..... CS, SC or degree vector
%    lmax .... Maximum degree of the Spherical Harmonic development
%
% OUT: 
%    vec ..... [l m Clm Slm] vector with Colombo ordering
%
% USES:
%    sc2cs

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2007-09-20: BD, initial version
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

[rows,cols] = size(mat);

if rows==lmax+1 && cols==rows 
    vec = (1:lmax+1)'*ones(1,lmax+1);
    ind = find(tril(vec)~=0);
    l   = vec(ind)-1;
    vec = vec';
    m   = vec(ind)-1;
    clm = mat(ind);
    mat(end,:) = [];
    mat = [zeros(1,lmax+1);mat]';
    slm = mat(ind);
    vec = [l m clm slm];
elseif rows==lmax+1 && 2*lmax+1==cols
    vec = (1:lmax+1)'*ones(1,lmax+1);
    ind = find(tril(vec)~=0);
    l   = vec(ind)-1;
    vec = vec';
    m   = vec(ind)-1;
    mat = sc2cs(mat);
    clm = mat(ind);
    mat(end,:) = [];
    mat = [zeros(1,lmax+1);mat]';
    slm = mat(ind);
    vec = [l m clm slm];
elseif rows==lmax+1 && cols==2
    vec = (mat(:,1)+1)*ones(1,lmax+1);
    ind = find(tril(vec)~=0);
    l   = vec(ind)-1;
    vec = vec';
    m   = vec(ind)-1;
    mat = (mat(:,2))*ones(1,lmax+1);
    clm = mat(ind);
    slm = clm;
    vec = [l m clm slm];
elseif rows==lmax+1 && cols==3
    vec = (mat(:,1)+1)*ones(1,lamx+1);
    ind = find(tril(vec)~=0);
    l   = vec(ind)-1;
    vec = vec';
    m   = vec(ind)-1;
    tmp = (mat(:,2))*ones(1,lmax+1);
    clm = tmp(ind);
    tmp = (mat(:,3))*ones(1,lmax+1);
    slm = tmp(ind);
    vec = [l m clm slm];
elseif rows==sum(1:lmax+1) && cols==4
    vec = mat;
    [vec(:,2) ind] = sort(vec(:,2));
    if isequal(mat(:,2),vec(:,2))
        display('Input already in Colombo ordering')
    end
    vec = [vec(ind,1) vec(:,2) vec(ind,3:end)];
else
    error('Input matrix is not in a compatible format (CS, SC, or [l Cl (Sl)])')
end
