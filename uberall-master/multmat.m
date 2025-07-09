function C = multmat(A, B, Atrans, Btrans)

% multmat mulitplicates two matrices, where each line of a matrix
% represents a rotation matrix at certain time point
%
% IN:
%    A ........ first matrix                                         [m, 9]  
%    B ........ second matrix                                        [m, 9] 
%    Atrans ... flag for using the transpose of A (default: false)   [bool] 
%    Btrans ... flag for using the transpose of B (default: false)   [bool] 
%
%               The multiplication A*B will be performed. This function is especially
%               designed to support the mulitplication of 3x3 rotation matirces.
%               Sometimes this must be done for a time series of rotation matrices.
%               Instead of using a for-loop it is possible to use this function. Of
%               course, since Matlab only supports calculations for 2D-matrices, the
%               format of the input matices must be altered.
%               The input format for the matrices A and B is as follows:
%                     A(:, 1) = A11(t)   for all t
%                     A(:, 2) = A12(t)   for all t
%                     A(:, 3) = A13(t)   for all t
%                     A(:, 4) = A21(t)   for all t
%                     A(:, 5) = A22(t)   for all t
%                     A(:, 6) = A23(t)   for all t
%                     A(:, 7) = A31(t)   for all t
%                     A(:, 8) = A32(t)   for all t
%                     A(:, 9) = A33(t)   for all t
%                     e.g. A(t) = [A11(t), A12(t), A13(t); A21(t), A22(t), ...
%                     Analog for B and C!
%
% OUT:
%    C ........ multiplication of A and B                            [m, 9]  

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias Weigelt (MW), DoGE, UofC                    
%    Matthias Roth (MR), GI, Uni Stuttgart               
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-05-15: MR, small help text brush-up
%    2013-02-07: MR, removed bug: 9x9 matrices would have been transposed
%                    automatically; assume that the matrix is correct in
%                    that case
%    2012-09-27: MR, added: possibility to transpose B
%    2003-03-25: MW, initial version
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

% constants:
narginchk(2, 4);

if ~exist('Atrans', 'var')
    Atrans = false;
end
if ~exist('Btrans', 'var')
    Btrans = false;
end

if ~islogical(Atrans), error('Atrans must be logical.'); end
if ~islogical(Btrans), error('Btrans must be logical.'); end

% input check
[mA, nA] = size(A);
[mB, nB] = size(B);

if (mA ~= mB) || (nA ~= nB)
    error('Matrices A and B have different dimensions!');
end

if (mA ~= 9) && (nA ~= 9)
    error('Wrong input format for matrix A');
elseif (mA == 9) && (nA ~= 9)
    warning('Matrix A is [9 x m] instead of [m x 9] --> I turn it for you.');
    A = A';
end

if (mB ~= 9) && (nB ~= 9)
    error('Wrong input format for matrix B');
elseif (mB == 9) && (nB ~= 9)
    warning('Matrix M is [9 x m] instead of [m x 9] --> I turn it for you.');
    B = B';
end

% sorting
if Atrans
    A1 = A(:, [1 4 7]);  % equal to A matrix 1. column
    A2 = A(:, [2 5 8]);  % equal to A matrix 2. column
    A3 = A(:, [3 6 9]);  % equal to A matrix 3. column
else
    A1 = A(:, 1:3);  % equal to A matrix 1. row
    A2 = A(:, 4:6);  % equal to A matrix 2. row
    A3 = A(:, 7:9);  % equal to A matrix 3. row
end

if Btrans
    B1 = B(:, 1:3);  % equal to B matrix 1. row
    B2 = B(:, 4:6);  % equal to B matrix 2. row
    B3 = B(:, 7:9);  % equal to B matrix 3. row
else
    B1 = B(:, [1 4 7]);  % equal to B matrix 1. column
    B2 = B(:, [2 5 8]);  % equal to B matrix 2. column
    B3 = B(:, [3 6 9]);  % equal to B matrix 3. column
end

% 1.row of C
C = zeros(size(A1,1),9);
C(:,1) = sum(A1 .* B1, 2);
C(:,2) = sum(A1 .* B2, 2);
C(:,3) = sum(A1 .* B3, 2);

% 2.row of C
C(:,4) = sum(A2 .* B1, 2);
C(:,5) = sum(A2 .* B2, 2);
C(:,6) = sum(A2 .* B3, 2);

% 3.row of C
C(:,7) = sum(A3 .* B1, 2);
C(:,8) = sum(A3 .* B2, 2);
C(:,9) = sum(A3 .* B3, 2);

