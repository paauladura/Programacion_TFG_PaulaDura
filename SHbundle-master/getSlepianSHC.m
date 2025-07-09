function [SHC3slepian,SL_INFO__]  = getSlepianSHC(GB_vectors, GB_values, GB_orders,GB_INFO__,select,alpha,beta,gamma)
% 
% GETSLEPIANSHC re-arranges the eigenvectors of Gruenbaum's eigenvalue 
% problem to the |S\C|-format of spherical harmonic coefficients and rotated 
% the coefficients -- if Euler angles are given -- to another location on the 
% sphere.
% 
% IN:   
%   GB_vectors: eigenvectors of Gruenbaum's eigenvalue problem 
%               = spherical harmonic coefficients of Slepian functions
%               (sparse matrix, special format)    
%   GB_values:  eigenvalues of Gruenbaum's eigenvalue problem 
%               = information about concentration within the cap 
%   GB_orders:  order of the spherical harmonics 
%               (<default: sorting with magnitude of GB_values!)
%   GB_INFO__:    cell array of further information: 
%                * .radius_rad ...size of spherical cap in radian  
%                * .lmax ........ maximum degree of expansion       
%                * .region....... shape of the region
%                * .algorithm.... of calculation: 'gruenbaum_evp'
%   select:.... index for the required Slepian basis functions      
%   alpha:..... angle of the first rotation around z-axis                  [rad][1x1]
%   beta:.....  angle of the second rotation around y-axis                 [rad][1x1]
%   gamma:....  angle of the third rotation around z-axis                  [rad][1x1]
%
% OUT:
%   SHC3slepian:spherical harmonic coefficicients in |C\S|-format and in a 3D matrix
%               (each layer of the 3D matrix is a individual Slepian basis function)
%                 
%   SL_INFO__:  cell array of further information: 
%                * .radius_rad ...size of spherical cap in radian  
%                * .lmax ........ maximum degree of expansion     
%                * .alpha_rad  .. angle of the 1. rotation around z-axis   [rad][1x1]
%                * .beta_rad  ... angle of the 2. rotation around y-axis   [rad][1x1]
%                * .gamma_rad ..  angle of the 3. rotation around z-axis   [rad][1x1]
%                * .eigenvalues ..Gruenbaum eigen values
%
%      
% USES: 
%    rotate_SHC
%
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
error(nargchk(4,8,nargin)); 
if nargin < 5,  select = 1;   end
if nargin < 6,  alpha = 0;    end
if nargin < 7,  beta = 0;     end
if nargin < 8,  gamma = 0;    end


if isnumeric(GB_vectors) == false
    error('Eigenvectors ''GB_vector'' must be a matrix of numerical values')
end 
if isnumeric(GB_values) == false
    error('Eigenvalues ''GB_values'' must be a vector of numerical values')
end  
if isnumeric(GB_orders) == false
    error('Order ''GB_orders'' must be a vector of numerical values')
end 
if isstruct(GB_INFO__) == false
    error('Information  ''GB_INFO__'' must be struct')
end 
if isnumeric(alpha) == false || isnumeric(beta) == false || isnumeric(gamma) == false
    error('all 3 Euler angles must be numerical values')
end
if isscalar(alpha) ==  false  || isscalar(beta) ==  false || isscalar(gamma) ==  false
    error('all 3 Euler angles must be scalar')
end  


%% preparation
lmax   = GB_INFO__.lmax;              % maximum degree of expansion
lmax1  = lmax+1;                      % auxilary value 
tmp    = fix(0.5*(lmax.^2+3*lmax))+1; % amount of Slepian functions in theory
select = intersect(0:tmp,select);     % restrict of 'select' to integer numbers
MAX    = numel(select)                % actual amount of Slepian functions
if MAX == 0,  error('Please insert meaningfull selection index');   end 
 
% initialization by zeros 
SHC3slepian(lmax1,lmax1,MAX) = 0;

for ii = 1:MAX
     % selection of spherical harmonic coefficients (spherical caps around Pole)
     m  = GB_orders(select(ii));
     cs = full(GB_vectors(m+1:lmax1,select(ii)));
     if m == 0                          % cosine coefficients
        SHC3slepian(:,1,ii) = cs;
     else                               % sine coefficients
        SHC3slepian(m,m+1:lmax1,ii) = cs;
     end
end
  

%% rotate to another location 
if (alpha^2+beta^2+gamma^2) > 0
    SHC3slepian  = rotate_shc(SHC3slepian,alpha,beta,gamma);
end

%% cell array with information about the Slepian basis functions
SL_INFO__             = GB_INFO__
SL_INFO__.alpha_rad   = alpha;  
SL_INFO__.beta_rad    = beta; 
SL_INFO__.gamma_rad   = gamma; 
SL_INFO__.eigenvalue  = GB_values(1,select); 
SL_INFO__.order       = GB_orders(1,select); 
 