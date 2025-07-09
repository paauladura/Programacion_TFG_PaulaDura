function C = multmatvek(A, B, Atrans)

% multmatvek mulitplicates a matrix and a vector, where each line of 
% them represents a rotation matrix at certain time point
%
% IN:
%    A .......... matrix                                                   [m, 9]  
%    B .......... vector                                                   [m, 3]  
%    Atrans ..... flag for using the transpose of A (default: false)       [bool]  
%
%                 The multiplication A*B will be performed. This function is especially
%                 designed to support the mulitplication of 3x3 rotation matrices.
%                 Sometimes this must be done for a time series of rotation matrices.
%                 Instead of using a for-loop it is possible to use this function. Of
%                 course, since Matlab only supports calculations for 2D-matrices, the
%                 format of the input matrices must be altered.
%                 The input format for the matrices A is as follows:
%                       A(:,1) = A11(t)   for all t
%                       A(:,2) = A12(t)   for all t
%                       A(:,3) = A13(t)   for all t
%                       A(:,4) = A21(t)   for all t
%                       A(:,5) = A22(t)   for all t
%                       A(:,6) = A23(t)   for all t
%                       A(:,7) = A31(t)   for all t
%                       A(:,8) = A32(t)   for all t
%                       A(:,9) = A33(t)   for all t
%                  e.g. A(t) = [A11(t), A12(t), A13(t); A21(t), A22(t), ...
%                  Analog for B!
%                       B(:,1) = B1(t)    for all t
%                       B(:,2) = B2(t)    for all t
%                       B(:,3) = B3(t)    for all t
%
% OUT:
%    C .......... multiplication of A and B                                 [m, 3]  

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias Weigelt (MW), DoGE, UofC
%    Matthias Roth (MR), GI, Uni Stuttgart     
%    Markus Antoni (MA), GI, Uni Stuttgart     
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-02-13: MA, consistent use of Atrans, whos for wrong dimensions
%    2013-02-08: MR, added error message if arrays A and B have different
%                    length, tiny help text brush up
%    2004-06-29: MW, initial version
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

narginchk(2, 3);

if ~exist('Atrans', 'var')
    Atrans = false;
end
if islogical(Atrans)== 0 || isscalar(Atrans) == 0
    error('Atrans must be logical and scalar.'); 
end

% input check
[rowA, columnA] = size(A);
[rowB, columnB] = size(B);



if (rowA ~= 9) && (columnA ~= 9)
    whos
    error('dimension:error','Wrong input format for matrix A');
elseif (rowA == 9) && (columnA ~= 9)
    warning('Array of matrices A is [9 x m] instead of [m x 9] --> I turn it for you.');
    A = A';
    rowA = columnA;
end

if (rowB ~= 3) && (columnB ~= 3)
    whos
    error('dimension:error','Wrong input format for matrix B');
elseif (rowB == 3) && (columnB ~= 3)
    warning('Array of vectors B is [3 x m] instead of [m x 3] --> I turn it for you.');
    B = B';
    rowB = columnB;
end


if (rowA ~= rowB) 
    whos
    error('dimension:error','Array of matrices A and array of vectors B have different length!');
end


% sorting
if Atrans
    disp('transpose matrix A')
    A1 = A(:, [1 4 7]);  % equal to A matrix 1. column
    A2 = A(:, [2 5 8]);  % equal to A matrix 2. column
    A3 = A(:, [3 6 9]);  % equal to A matrix 3. column
else
    A1 = A(:, 1:3);  % equal to A matrix 1. row
    A2 = A(:, 4:6);  % equal to A matrix 2. row
    A3 = A(:, 7:9);  % equal to A matrix 3. row
end

% 1.row of C
C1 = sum(A1 .* B, 2);
C2 = sum(A2 .* B, 2);
C3 = sum(A3 .* B, 2);

% form new format
C = [C1 C2 C3];
