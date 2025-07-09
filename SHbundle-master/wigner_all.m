function [dlmkBlock,FxlmkBlock, order_m,order_k] = wigner_all(lmax,inc,output)
% [dlmkBlock,FxlmkBlock, order_m,order_k] = wigner_all(lmax,inc,output)
% 
% WIGNER_ALL caluclates either 
%   * all Wigner-d-functions 'dlmk' of the irreducible representation 
%     of the rotation groups SU(2) and SO(3) up to degree l = 'lmax' 
%   * or all inclination function 'flmk' and 'fxlmk' of satellite geodesy
%     up to degree 'lmax', i.e. the  product of Wigner-d-functions and the 
%     Legendre functions at the equator
%
% Based on the Wigner-d-functions, the complex inclination functions are 
% calculated via
%
%   flmk(inc) = i^(m-k) *  dlmk(inc) * P_{l,m}(0)
%   fxlmk(inc) = i^(m-k) *  dlmk(inc) * ( d{P_{l,m}(0)} / d{theta} )
%
% where P_{l,m}(0) are the Legendre functions (at the equator) in complex 
% normalization. 
%
% IN:
% lmax:...... maximum degree of the Wigner-d-functions                          [1x1]
% inc:....... angle of rotation (scalar)                                   [rad][1x1]
% output:.... change output variable                                            [string]
%       *   'flmk':    first output is inclination function in complex form
%                      (2nd output: crosstrack inclination functions Fxlmk )
%       *   'dlmk':    first output is dlmk function of SO(3) group
%                      (2nd output: empty matrix)       
%
% OUT
% dlmkBlock:... Wigner-d-function 'dlmk' or inclination function 'flmk'         [cell quadratic matrices]   
%               arrangement of the dlmk/flmk 'in a pyramide'
%                   - cell-index/level j refers to degree j = l+1
%                   - each matrix in level j contains all function of order -l <=(m,k)<=l    
%                   - sorting in the matrix accoding to order_k/order_m output      
% FxlmkBlock:.. crosstrack inclination functions Fxlmk or empty                 [cell quadratic matrices]   
%               arrangement of the fxlmk 'in a pyramide'
%                   - cell-index/level j refers to degree j = l+1
%                   - each matrix in level j contains all function of order -l <=(m,k)<=l 
%                   - orting in the matrix accoding to order_k/order_m output 
% order_m:..... 1. order of the last dlmk (in the cell{lmax+1})
% order_k:..... 2. order of the last dlmk (in the cell{lmax+1})
%     
% REMARKS:
%    *   There are also algorithms, which calculate the inclination
%        function in one recursion without final multiplication dlmk * plm
% 
% LITERATURE:
%        P. J. KOSTELEC & D. N. ROCKMORE (2003): 
%        "FFTs on the Rotation Group"
%        P. J. KOSTELEC & D. N. ROCKMORE (2003): 
%        "SOFT: SO(3) Fourier Transforms"
%
% USES:
%    Legendre0

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2021-04-12: MA, brushed up
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

error(nargchk(2,3,nargin)); 

if nargin < 3
    output = 'dlmk';
end

%% Check of input arguments:
if isscalar(inc)~=1 
    error('only scalar rotation angles ''inc'' are possible')
end
if max(abs(inc))>pi
    warning('rotation angle ''inc'' must be given in radian')
end;

if isscalar(lmax)~=1 
    error('maximum degree ''lmax'' must be scalar')
end
if lmax<0 || rem(lmax,1)~=0
    error('maximum degree must be a positive integer')
end
if isempty(output)== true
    output = 'dlmk';
end

%% allocation of memory
for degree_l = 0:lmax
    dlmkBlock{degree_l+1} = eye(degree_l*2+1);
end

if inc == 0
    warning('the angle ''inc = 0'' produces only identity matrices for dlmk; NO CALCULATION of flmk and fxlmk')
    FxlmkBlock = [];
    return
end

%% dimensions 
mmax = lmax+1;
nmax = 2*lmax+1;
nmax1= nmax+1;



%% allocation of 2nd output
FxlmkBlock = dlmkBlock;

Dm2(nmax,nmax)    = 0;
Dm1(nmax,nmax)    = 0;
Wnkm(nmax,1)      = 0;
powerCosZ(nmax,1) = 0;
powerSinZ(nmax,1) = 0;


%% trigonometric functions
cosb      = cos(inc);
powerCosA = cos(inc/2).^(nmax:-1:0); 
powerSinA = sin(inc/2).^(0:nmax); 


%% order and factorials (logarithm calculation of the norm)
[order_m,order_k] = meshgrid([-lmax:lmax]);
km     = order_m.*order_k;
mm_sqr = order_m.^2;
kk_sqr = order_k.^2;
kv     = order_k(:,1);
LOGJ(mmax:nmax+lmax) = (gammaln(2:nmax+1))/2;



%% Legendre functions
switch output
case 'flmk'
     fprintf(1,'\n WIGNER_ALL.M: calculation of the (complex) incliation functions ''flmk'' and ''fxlmk''')
    [Pnm0e,dPnm0e] = Legendre0(lmax,'cplusminus');
    Pnm0e   = Pnm0e';
    dPnm0e  = dPnm0e';
    ikm     = 1i.^(order_m-order_k);
case 'dlmk'
    fprintf(1,'\n WIGNER_ALL.M: calculation of the Wigner-d-functions ''dlmk''')
    Pnm0e   =  ones(1,nmax);
    dPnm0e  = ones(1,nmax);
    ikm     = 1;
otherwise
    error('sorry this option is not implemented')
end   


 
%% index and first values 
index1 = mmax;
index2 = mmax;

dlmkBlock{1} = 1;
FxlmkBlock{1} = 0;
Dm1(index1,index1) = 1;
signum = -(-1).^kv;
plus1 = 1;
count = 3;

for degree_l = 0:lmax-1
    
    %% 3-term recursion (unnormalized according to SOFT [17])
    degree_plus = degree_l+1;
    degree_sqr  = degree_l^2;
    S1 = sqrt(abs(degree_plus^2-mm_sqr));
    S2 = sqrt(abs(degree_plus^2-kk_sqr));
    a  = degree_plus./S1 .* (2*degree_l+1)./S2 .*Dm1;
 
    b  = sqrt(abs(degree_sqr - mm_sqr))./S1 .* sqrt(abs(degree_sqr - kk_sqr))./S2 .* Dm2;
    Dm0= a .* (cosb-km./((degree_l+plus1)*degree_plus)) - b *(degree_plus)/(degree_l+plus1);
    
   
    % update of index
    index1 = index1-1;
    index2 = index2+1;
    


    %% correct margin of the squares:
    % factor = sqrt ( (2*J)! / [(J+M)! * (J-M)!] ) in log-calculation
    Wnkm(:,1) = exp((LOGJ(mmax+1+2*degree_l) - LOGJ(lmax+degree_plus+kv) - LOGJ(lmax+degree_plus-kv)));

    powerSinZ(index1:index2) = powerSinA(1:count);
    powerCosZ(index1:index2) = powerCosA(nmax-count+2:nmax1);
    Dm0(:,index1) = Wnkm.*powerCosZ.*powerSinZ;    % dJ[-J,M]
    Dm0(index1,:) = Dm0(:,index1).*signum;         % dJ[M,-J]
    Dm0(index2,:) = flipud(Dm0(:,index1));         % dJ[M,J]
    Dm0(:,index2) = Dm0(index2,:).*signum';        % dJ[J,M]
  
    %% multiply with Legendre functions to get inclination (or with 1 to keep dlmk)
    fnmk   = ikm .* (Dm0*diag(Pnm0e(:,degree_plus+1)));
    fxnmk  = ikm .* (Dm0*diag(dPnm0e(:,degree_plus+1)));
    
    %% store into cell array 
    dlmkBlock{degree_plus+1}  = fnmk(index1:index2,index1:index2);
    FxlmkBlock{degree_plus+1} = fxnmk(index1:index2,index1:index2); 
    
    %% updates
    count = count+2;  
    Dm2 = Dm1;
    Dm1 = Dm0;
    signum = -signum;
    plus1 = 0;
    
end

switch output
case 'dlmk'
    FxlmkBlock = [];
otherwise
end

end

