function [index] = sortLegendre(kmax,FromTo)
% SORTLEGENDRE provides an index for re-arranging the Legendrefunctions 
% their derivatives, the spherical harmonic coefficents or other related 
% functions between different 'sorting styles'.
%
% The standard sortings are the DEGREE-ORDER format or the ORDER-DEGREE
% format, where the functions are sorted in the following way:
%
% DEGREE-ORDER:
%           | f{0,0}(x1)   f{0,0}(x2)   f{0,0}(x3) ....   f{0,0}(xT)|
%           | f{1,0}(x1)   f{1,0}(x2)   f{1,0}(x3) ....   f{1,0}(xT)|
%           | f{1,1}(x1)   f{1,1}(x2)   f{1,1}(x3) ....   f{1,1}(xT)|
%           | f{2,0}(x1)   f{2,0}(x2)   f{2,0}(x3) ....   f{2,0}(xT)|
% F   =     | f{2,1}(x1)   f{2,1}(x2)   f{2,1}(x3) ....   f{2,1}(xT)|
%           |     :                                             :   |
%           | f{n,m}(x1)             ..                         :   |
%           |     :                                             :   |
%           | f{N,N-1}(x1) f{N,N-1}(x2) f{N,N-1}(x3) .... f{N,N}(xT)|
%           | f{N,N}(x1)   f{N,N}(x2)   f{N,N}(x3) ....   f{N,N}(xT)|
% ORDER-DEGREE:
%           | f{0,0}(x1)   f{0,0}(x2)   f{0,0}(x3) ....   f{0,0}(xT)|
%           | f{1,0}(x1)   f{1,0}(x2)   f{1,0}(x3) ....   f{1,0}(xT)|
%           | f{2,0}(x1)   f{2,0}(x2)   f{2,0}(x3) ....   f{2,0}(xT)|
% F   =     |     :                                             :   |
%           | f{n,0}(x1)             ..                         :   |
%           | f{1,1}(x1)   f{1,1}(x2)   f{1,1}(x3) ....   f{1,1}(xT)|
%           | f{2,1}(x1)   f{2,1}(x2)   f{2,1}(x3) ....   f{2,1}(xT)|
%           |     :                                             :   |
%           | f{n,1}(x1)             ..                         :   |
%           | f{2,2}(x1)   f{2,2}(x2)   f{2,2}(x3) ....   f{2,2}(xT)|
%           |     :                                             :   |
%           | f{N,N}(x1)   f{N,N}(x2)   f{N,N}(x3) ....   f{N,N}(xT)|
% where
%   f{n,m}(t) is the function of degree n, order m at location t.
%
% IN:
%    kmax ....... maximum degree of the Legendre functions                   [1x1]
%    FromTo ..... conversion between ...                                     ['sring'] 
%           * 'degree2order': DEGREE-ORDER -> ORDER-DEGREE
%           * 'order2degree': ORDER-DEGREE -> DEGREE-ORDER 
%
%
% OUT:
%    index ...... vector of re-arranging
%
% EXAMPLE:
%    kmax = 10
%    theRAD = 0:.1:pi;
%    % Legendrefunctions in DEGREE-ORDER format
%    [Pnm_DO, delVar,delVar,DegreeOrder]= legendreP(kmax,theRAD);
%    [indexDO] = sortLegendre(kmax,'degree2order')
%    % re-arrangement in ORDER-DEGREE
%    OrderDegree = DegreeOrder(indexDO,:)
%    Pnm_OD = Pnm_DO(indexDO,:);
%    % back to DEGREE-ORDER    
%    [indexOD] = sortLegendre(kmax,'order2degree');
%    DegreeOrder2 = OrderDegree(indexOD,:)
%    Pnm_DO2 = Pnm_OD(indexOD,:);
%
% SEE ALSO:
%    legendreIndex
%
% USES:
%    --

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% author:               
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-10-28: MA, initial version
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

%% check the degree argument
if isscalar(kmax) == 0
    error('maximal degree ''kmax'' must be scalar')
end
if kmax<0 || rem(kmax,1) ~= 0
    error('maximal degree ''kmax'' must be positive and integer')
end

%% initialization
num = (kmax.*(kmax+1))./2+(kmax+1); 
index(1:num,1) = 1;


%% index determination
switch lower(FromTo)
case 'degree2order'
    %% DEGREE-ORDER -> ORDER-DEGREE
    update = cumsum(0:kmax+1);
    in1    = 1; 
    last   = 1; 
    for ii = 1:kmax+1;
        in2     = in1+kmax-ii+1;
        index(in1:in2,1) = last+update(ii:kmax+1);
        last    = last+1;
        in1     = in2+1;    
    end
case 'order2degree'
    %% ORDER-DEGREE -> DEGREE-ORDER 
    update = cumsum([kmax:(-1):0]);
    in1    = 2;
    last   = 2;
    for ii = 1:kmax;
        in2     = in1+ii;
        index(in1:in2) = [last, last+update(1:ii)];
        last    = last+1;
        in1     = in2+1;
    end
otherwise 
    error(['sorting ''' FromTo '''is not implemented'])
end




