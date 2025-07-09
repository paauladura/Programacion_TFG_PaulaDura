function d = dlmk(l,beta)

% DLMK(L,BETA) creates a d-matrix for degree l and rotation angle beta.
%    The matrix has size (2l+1)*(2l+1). It is an irreducible matrix
%    representation of weight (=label=degree) l of a rotation about the y-axis.
%
%    INPUT: l    - spherical harmonic degree (weight or label)
%           beta - rotation angle (in radians)
%
%    First the matrix representation D(Ly) of an infitesimal rotation (Ly) 
%    is created. Ly is related to the step operator J+ and J- by:
%         Ly = -i Jy = -i( -i/2(J+ - J-) ) = 1/2(J- - J+) 
%    The d-matrix is expm(beta*Ly) then.
%
% Note: The orders start +l --> -l rather than the normal convention of
% -l --> +l.
%
% Disclaimer: It is not guaranteed that (i) the proper minus signs are there.
% Perhaps a (-1)^(m-k) factor should be used, and (ii) that the routines
% perform correctly for higher degrees. 
%

% Nico Sneeuw, Munich 17/03/94

dly  = zeros(2*l+1);
for i = 1:2*l
   m = l-i+1;
   dly(i+1,i) = sqrt((l+m)*(l-m+1));
end
dly = (dly - dly')/2;

d = expm(beta*dly); 

