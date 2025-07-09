% uberall toolbox - commonly used tools, v.2014
%
%   checkcoor       - preperation of spherical coordinates (unit, dimension)
%   diagindx        - provides the indices of the diagonal elements of a 
%                     given square matrix
%   getopt          - optional parameter handling by string tags
%   grule           - gauss points and weights for one-dimensional quadrature
%   isint           - test for being integer
%   ispec           - IFFT for real-valued Fourier spectrum
%   ispec2          - 2-D IFFT for real-valued Fourier spectrum
%   jump2NaN        - inserts NaN's when differences within a vector
%                     exceed a certain threshold (in absolute sense)
%                     (e.g. switch in longitude at +/-180 degree, or short arc
%                     analysis)
%   long360toeast   - converts longitude values from 0 to 360 to -180 to 180
%   lying           - put matrix "horizontally"
%   multmatmult     - mulitplicates multiple matrices in [9 x N] format 
%                     from the left to right
%   multmatvek      - multiplication of a set of vectors and matrices (in the
%                     terms of a [Nx9] respectively [Nx3] hypermatrix)
%   multmat         - multiplication of 2 rotation matrices in the [Nx3] format 
%   normg           - normal gravity (formula of Somigliana)
%   Q1              - 3D mirroring at X-axis
%   Q2              - 3D mirroring at Y-axis
%   Q3              - 3D mirroring at X-axis
%   R1              - 3D rotation around X-axis
%   R2              - 3D rotation around Y-axis
%   R3              - 3D rotation around Z-axis
%   replace         - replace elements in vector/matrix
%   spec            - FFT for real-valued Fourier spectrum
%   spec2           - 2-D FFT for real-valued Fourier spectrum
%   standing        - put matrix "upright" 
%   trapstrip       - strip a matrix off its corners
%   twaitbar        - displays a text waitbar in the command window
%   xyz2mat         - converts a look-up table for matrix elements to matrix
%
%   constants       - Constants like GM, semi-major axis of the earth
%   constants_grace - Constants specific to GRACE satellite
%   constants_champ - Constants specific to CHAMP satellite
% 
% Data:
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% authors: 
%   Nico Sneeuw (NS)
%   Matthias Weigelt (MW)
%   Markus Antoni (MA)
%   Balaji Devaraju (BD)
%   Matthias Roth (MR)
%
% file compilation: 
%   v.2014:  2014-10-23, MA, last update of 2014
%   v.2013:  2013-02-13, MR, moving the general tools from SHbundle to new "uberall"
%                            toolbox

