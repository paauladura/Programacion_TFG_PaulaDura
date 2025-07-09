function SHC = rotate_shc(shc,alpha,beta,gamma)
% 
% Rotation of the spherical harmonic coefficients, so that the a rotated field 
% quantity can be obtained by a summation over the original SHS and
% the transformed coefficients.
% 
% shc: ..... spherical harmonic coefficients in |C\S|-format
%            either in real or complex geodetic normalization 
%            (multiple sets of coefficients can be given in 3D matrix)
% lmax: .... maximum degree of the expansion
% alpha:.... angle of the first rotation around z-axis                     [rad][1x1]
% beta:..... angle of the second rotation around y-axis                    [rad][1x1]
% gamma:...  angle of the third rotation around z-axis                     [rad][1x1]
%
% SHC: .....  spherical harmonic coefficients 
%               (same format/normalization as the input argument)
%
% USES:
%    real2cpxsh, wigner_all, cpx2realsh, cs2sc, sc2cs

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2021-04-09: MA, brushed up
%    2021-04-01: MA, initial version
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

error(nargchk(4,4,nargin)); 

%% Checking
% ... for dimension
if max([numel(alpha),numel(beta),numel(gamma)])>1
    error('only scalar angles ''alpha'', ''beta'', ''gamma'' are possible')
end
% ... for unit/type
if max([alpha,beta,gamma]) > 2*pi
  warning('Are all angles given in radian?')
end


  
%% Checking spherical harmonic coeffcients
%% (sometimes it can be helpful to apply the same rotation to several sets of coefficients)
[lmax1,col,depth] = size(shc);
lmax = lmax1-1;
if col == lmax1
    if col < 7
       warning('assuming the |C\S|-format; If not correct, increase the matrix by zero padding')
    end
else
    error('please insert SH coefficients in |C\S|-format')
end

%% all-Wigner-functions by recursive computation
dlmkBlock = wigner_all(lmax,beta,'dlmk');

%% initialization 
SHC(lmax1,col,depth) = 0; 
exp_gamma = (exp(-1i*[-lmax:lmax]*gamma));
exp_alpha = (exp(-1i*[-lmax:lmax]*alpha));
for ii = 1:depth 
    shc_i = shc(:,:,ii);
    if isreal(shc_i) == true  
        fprintf(1,'\n ROTATE_SHC.M: converting coefficients from real to complex representation\n');
        cis = real2cpxsh(shc_i,lmax);
        BoolianComplex = false
    else
        BoolianComplex = true
        cis = shc_i;
    end

    % rotation according to SO(3) group
    cis_rot =  bsxfun(@times, cs2sc(cis), exp_gamma);
    for ll = 0:lmax
        d = dlmkBlock{ll+1};
        cis_rot(ll+1,lmax1-ll:lmax1+ll) = cis_rot(ll+1,lmax1-ll:lmax1+ll)*d;
    end
    cis_rot =  bsxfun(@times,cis_rot,exp_alpha)  ;
    %% back to real SH-coeffcients if necessary
    if BoolianComplex == false;
        fprintf(1,'\n ROTATE_SHC.M: inverse converting coefficients from complex to real representation\n');
        cis_rot = cpx2realsh(cis_rot ,lmax);
    end
    cis_rot = sc2cs(cis_rot);
    SHC(:,:,ii) = cis_rot;
end

