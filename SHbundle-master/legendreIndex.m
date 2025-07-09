function index = legendreIndex(n,m,format)

% LEGENDREINDEX determines the rows to identify certain functions in the 
%       [1] "degree-order format" (default) 
%       [2] or in the "order-degree format". 
%
% Add [1]:
% Several programs (e.g. legendreP.m, wigner.m,...) provide a matrix of
% all functions up to degree N. The output is arranged in a matrix
%
%           | f{0,0}(x1)   f{0,0}(x2)   f{0,0}(x1) ....   f{0,0}(xT)|
%           | f{1,0}(x1)   f{1,0}(x2)   f{1,0}(x1) ....   f{1,0}(xT)|
%           | f{1,1}(x1)   f{1,1}(x2)   f{1,1}(x1) ....   f{1,1}(xT)|
%           | f{2,0}(x1)   f{2,0}(x2)   f{2,0}(x1) ....   f{2,0}(xT)|
% F   =     | f{2,1}(x1)   f{2,1}(x2)   f{2,1}(x1) ....   f{2,1}(xT)|
%           |     :                                             :   |
%           | f{n,m}(x1)             ..                         :   |
%           |     :                                             :   |
%           | f{N,N-1}(x1) f{N,N-1}(x2) f{N,N-1}(x1) .... f{N,N}(xT)|
%           | f{N,N}(x1)   f{N,N}(x2)   f{N,N}(x1) ....   f{N,N}(xT)|,
%
% where the first number in the {}-brackets denotes the degree n and the
% second the order m (with 0<=m<=n<=N). In some cases only the functions of
% a certain degree or certain order is necessary. 
%
%    example: 
%    DegreeOrder(:,1) = [0 1,1 2,2,2  3,3,3,3  4,4,4,4,4  5,5,5,5,5,5];
%    DegreeOrder(:,2) = [0 0,1 0,1,2  0,1,2,3  0,1,2,3,4  0,1,2,3,4,5];
%    index_1  = legendreIndex(5,5)
%    degree5_5 = DegreeOrder(index_1,:) 
%    degree4_123 = DegreeOrder(legendreIndex(4,[1,2,3]),:) 
%    degree345_2 = DegreeOrder(legendreIndex([3,4,5],2),:) 
%
% Add [2]:
% Speed optimized functions like Legendre_mex provide the sorting according
% to order and then to degree
%           | f{0,0}(x1)   f{0,0}(x2)   f{0,0}(x1) ....   f{0,0}(xT)|
%           | f{1,0}(x1)   f{1,0}(x2)   f{1,0}(x1) ....   f{1,0}(xT)|
%           | f{2,0}(x1)   f{2,0}(x2)   f{2,0}(x1) ....   f{2,0}(xT)|
% F   =     |     :                                             :   |
%           | f{n,0}(x1)             ..                         :   |
%           | f{1,1}(x1)   f{1,1}(x2)   f{1,1}(x1) ....   f{1,1}(xT)|
%           | f{2,1}(x1)   f{2,1}(x2)   f{2,1}(x1) ....   f{2,1}(xT)|
%           |     :                                             :   |
%           | f{n,1}(x1)             ..                         :   |
%           | f{2,2}(x1)   f{2,2}(x2)   f{2,2}(x1) ....   f{2,2}(xT)|
%           |     :                                             :   |
%           | f{N,N}(x1)   f{N,N}(x2)   f{N,N}(x1) ....   f{N,N}(xT)|,
%
% in this case, the 1. argument is interpreted as lmax and the 2. argument 
% as the desired (scalar) degree
%
%    example: 
%    OrderDegree(:,1) = [0,1,2,3,4,5  1,2,3,4,5  2,3,4,5  3,4,5 4,5  5];
%    OrderDegree(:,2) = [0,0,0,0,0,0  1,1,1,1,1  2,2,2,2  3,3,3 4,4  5] 
%    index_2  = legendreIndex(5,5,'orderdegree')
%    degree5 = OrderDegree(index_2,:) 
%    degree2 = OrderDegree(legendreIndex(5,2,'orderdegree'),:) 
%    degree1 = OrderDegree(legendreIndex(5,1,'orderdegree'),:) 
%    degree0 = OrderDegree(legendreIndex(5,0,'orderdegree'),:) 
%
% IN:
%    n ........ [Nx1]    degree of the function                   
%    m ........ [Nx1]    order of the function [1]/ maximal degree [2]                          
%    format..[string]    sorting method
%                        'degreeorder' (default)
%                        'orderdegree'
% OUT:  
%    index .... [Nx1]    rows of the degree-order/order-degree format              
%
% REMARKS:
%    Inb the degree-order format, one of the inputs can be a vector, while 
%    the second must be scalar!
%    the maximal order m is smaller or equal the degree n
%    in geodesy degree and order must be positive integers, and so 
%    possible fractional or negative orders lead to an error in the
%    latter use

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-09-26: MA, added order-degree format
%    2013-01-21: MA, last modification
%    2011-03-23: MA, initial version
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

if nargin<3
    format = 'degreeorder' ;
end

n(rem(n,1)~=0) = [];
m(rem(m,1)~=0) = [];
if nargin<2
    if min(n)<0
        error('nonInteger:error','degree must be non negative')
    end
    index = fix(0.5*(n.^2+3*n))+1;
    return
end

if min(n)<0 || min(m)<0
    error('nonInteger:error','degree and order must be non negative')
end



switch lower(format)
case 'degreeorder' 

    ii = m>max(n);
    m(ii) = [];
    ii = n<min(m);
    n(ii)= [];

    if min(length(n),length(m))~=1
        error('nonInteger:error','degree or order must be scalar')
    end
    index = (n.*(n+1))./2+(m+1); 
case 'orderdegree'
    lmax = n;
    ii = m>lmax;
    m(ii) = [];
    if length(n)+length(m)~=2
        error('nonInteger:error','degree and maximal degree must be scalar')
    end
    index = (m+1)+ [0, cumsum(lmax:(-1):(lmax-m+1))]; 
otherwise
    error('sorry this sorting is not implemented yet')
end
