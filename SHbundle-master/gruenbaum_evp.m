function [GB_vectors, GB_values, GB_orders, GB_INFO__]= gruenbaum_evp(Theta0,lmax,SortEigen)
% [GB_vectors, GB_values, GB_orders, GB_INFO__]= gruenbaum_evp(Theta0,lmax,SortEigen)
%
% GRUENBAUM_EVP solves the eigenvalue problem of the Gruenbaum matrix [T], i.e.
%
%       [T] * gm> = chi * gm>
%
% where order-dependent eigen vector gm> provides the spherical harmonic 
% coefficients of the Slepian basis functions within a SPHERICAL CAP. 
%
% Please note, that the standard caluclation 
%
%      [D] * gm> = Lambda * gm>
%
% has the same eigenvectors (i.e. coefficients), but different eigenvalues. 
%  
% IN:
%   Theta0:.... co-latitude (i.e. spherical distance) to the North-pole    [rad][1x1] 
%   lmax:...... maximum degree of expansion                                 [--][1x1] 
%   SortEigen:re-ordering of the eigenvectors                                  ['string']
%           *  sorted w.r.t. to magnitude of eigenvalues <default>
%           *  original sorting per order m  (here: SortEigen = 'keepsorting'>)
% 
% OUT:
%   GB_vectors: eigenvectors of the Gruenbaum eigenvalue problem 
%               = spherical harmonic coefficients of Slepian basis functions
%               (sparse matrix, special format)    
%   GB_values:  eigenvalues of the Gruenbaum eigenvalue problem 
%               = information about concentration within the cap 
%   GB_orders:  order of the spherical harmonics 
%               (<default: sorting with magnitude of GB_values!)
%   GB_INFO__:    cell array of further information: 
%                * .radius_rad ...size of spherical cap in radian  
%                * .lmax ........ maximum degree of expansion       
%                * .region....... shape of the region
%                * .algorithm.... of calculation: 'gruenbaum_evp'
% 
% ATTENTION
%      *    The eigenvalue problem is solved per order m. For memory and understanding,  
%           the coefficients have a special format in a sparse matrix. Please apply the 
%           routine  < GETSLEPIANSHC > for a standard |C\S|-format per Slepian basis
%           function.
%
% LITERATUR: 
%       *   Simons & Dahlen: Spherical Slepian functions and the polar gap in geodesy
%           Geophys. J. Int 2006, 1039-1061

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2021-10-20: MA, initial version
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


%% checking the input 
error(nargchk(2,3,nargin)); 

% sorting of eigenvalues and eigenvectors
if nargin < 3,  SortEigen = 'unchanged';  end
if ischar(SortEigen) == false
    error('re-arrangement ''SortEigen'' requires string-input')
end
% checking the co-latitude <<Theta0>>
if isnumeric(Theta0) == false
     error('co-latitude ''Theta0'' must be numerical value ')
end 
if isscalar(Theta0) == false 
    error('co-latitude ''Theta0'' must be scalar')
end
if Theta0<0 || Theta0>pi
    error('co-latitude ''Theta0'' must be given in radian [0,pi])')
end
% checking the maximum degree lmax
if isnumeric(lmax) == false
    error('maximum degree '' lmax'' must be a numerical value ')
end
if isscalar(lmax) == false 
    error('maximum degree ''lmax'' must be scalar')
end
if rem(lmax,1) ~= 0
    error('maximum degree ''lmax'' must be integer')
end



%% auxilary values
lmax1 = lmax+1;
ll    = 0:lmax;

%% allocation
MAX   = fix(0.5*(lmax.^2+3*lmax))+1;
GB_vectors  = spalloc(lmax1,MAX, sum((ll+1).^2));
GB_values(1,MAX) = NaN;
GB_orders(1,MAX) = NaN;

%% Gruenbaum's eigenvalue problem
% elements on the diagonal (formula 79a)
D = diag([-ll.*(ll+1).*cos(Theta0)]);
ind1 = 1;
for mm = 0:lmax   
    % parallels to diagonal (formula 79b)
    llm  = ll(mm+1:lmax);
    t_pm = [llm.*(llm+2) - lmax*(lmax+2)] .* sqrt([(llm+1).^2-mm^2]./[(2*llm+1).*(2*llm+3)]);
    T    = D(mm+1:lmax1,mm+1:lmax1) + diag(t_pm,1) +diag(t_pm,-1);

    % EVP per matrix 
    [Gm,chi]=eig(T);
  
    % indices and saving
    ind2  = ind1 +lmax-mm;
    GB_vectors(mm+1:lmax1,ind1 :ind2) = Gm;
    GB_values(1,ind1 :ind2) = diag(chi);
    GB_orders(1,ind1 :ind2) = mm;
    ind1  = ind2 +1;
end

%% sorting of the eigenvalues and eigenvectors
switch lower(SortEigen)
case 'keepsorting'
    disp('keeping the original sorting of GRUENBAUM eigenvalues')
otherwise
    disp('sorting of eigen vectors according to magnitude of GRUENBAUM eigenvalues')
    [tmp, index] = sort(GB_values(1,:));
    GB_values  = GB_values(index);
    GB_orders  = GB_orders(index);
    GB_vectors = GB_vectors(:,index);  
end
    
%% cell array with information about the eigenvalue problem
GB_INFO__.radius_rad = Theta0;    
GB_INFO__.lmax = lmax;
GB_INFO__.region    = 'spherical cap';
GB_INFO__.algorithm = 'gruenbaum_evp';
