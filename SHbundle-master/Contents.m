% Spherical Harmonic Computation and Graphics tools, v.2021
% Requires uberall toolbox
%
%-----------------------------------------------------------------------
% SH Synthesis & Analysis:
%-----------------------------------------------------------------------
%   eigengrav     - isotropic spectral transfer (eigenvalues)
%   gradpshs      - Global Spherical Harmonic Synthesis - gradient
%   gsha          - Global Spherical Harmonic Analysis
%   gshs_         - Global Spherical Harmonic Synthesis
%   gshs_grid     - Global Spherical Harmonic Synthesis - on any grid
%   gshs_ptw      - Global Spherical Harmonic Synthesis - pointwise
%   lovenr        - Love-numbers for the deformation of the elastic Earth
%   nablaPot      - radial derivatives for the disturbing potential along
%                   the orbit 
%   pshs          - pointwise Spherical Harmonic Synthesis
%   regionalshs   - Regional spherical harmonic synthesis on regular grids
%   shspkamph     - Amplitude and phase of spherical harmonic spectrum
%   sph_gradshs   - creates a SINGLE POINT HANDLE for calculating the
%                   gradient of a potential
%   tenspshs      - pointwise Spherical Harmonic Synthesis of the Marussi
%                   tensor 
%   upwcon        - upward continuation in the spectral domain
%   updwcont      - pointwise upward continuation in the spatial domain
%-----------------------------------------------------------------------
% SH Synthesis of two-point functions on the sphere:
%-----------------------------------------------------------------------
%   covord        - converts covariance matrix into a four variable cell
%                   structure
%   gshs2ptfun    - Same as GSHSCOV but used for computing the covariance
%                   function/bi-polar field for a particular point on the
%                   sphere iand for any type of grid (also individual points).
%   gshscov       - Synthesizes SH spectral covariance matrices to spatial
%                   covariances
%-----------------------------------------------------------------------
% Legendre polynomials & SH:
%-----------------------------------------------------------------------
%   diffLegpol     - unnormalized Legendre polynomials and their 1., 2. and
%                   3. derivatives  
%   iplm          - integrated Legendre functions
%   Legendre0     - normalized Legendre functions at the equator, with
%                   exact zeros of odd functions and in the triangle format
%                   (for construction the inclination functions via products)
%   Legendre_mex  - fully normalized Legendre functions of first kind up to
%                   degree n (order of coefficients different to 
%                   legendreP); uses X-numbers for a stable computation
%                   from degree/order >= 1500.
%   legendreP     - all Legendre functions up to degree n and the 1. 
%                   and 2. derivatives 
%   legendreP_Xnum - the same as legendreP, recommended for degree/order > 1500
%   legpol        - unnormalized Legendre functions
%   neumann       - weights and nodes for Neumann's quadrature methods
%   plm           - Legendre functions (normalized) and derivatives, consider
%                   using Legendre_mex for speed!
%   ylm           - surface spherical harmonics
%-----------------------------------------------------------------------
% Visualization:
%-----------------------------------------------------------------------
%   shplot        - plot SH spectrum
%   shprepare     - prepare SH spectrum for plotting purposes
%   ylmplot       - plot surface spherical harmonic
%-----------------------------------------------------------------------
% SH coefficients storage format conversion:
%-----------------------------------------------------------------------
%   checkshformat - Checks the format of the input SH coefficients
%   clm2klm       - Convert order-wise look-up-table format to degree-wise
%   clm2sc        - Indexed column-vector-format to SC-format, with 
%                   maximum degree and order
%   cs2sc         - CS-format to SC-format conversion (consult
%                   documentation for format explanations) 
%   cs2vec        - CS-format to column-vector-format
%   cssc2clm      - CS/SC-format to indexed-column-vector format
%   klm2clm       - Converts back from degree-wise look-up-table format
%                   to order-wise
%   legendreIndex - indices to select spherical harmonic coefficients of
%                   certain degree/order 
%   parse_icgem   - Reads SH coefficients from ICGEM-format ASCII files
%                   and returns them in indexed column-vector-format
%   sc2cs         - SC-format to CS-format conversion
%   sortLegendre  - conversion between degree-order format and 
%                   order-degree format.
%   vec2cs        - Column-vector-format to CS-format%
%-----------------------------------------------------------------------
% Statistics:
%-----------------------------------------------------------------------
%   degcorr       - SH degree correlation  
%   degreeinfo    - SH information per degree (signal or error degree RMS;
%                   commission error) 
%   globalmean    - Computes the global mean of the field
%   globalpower   - Computer the power of a given field
%   kaula         - Kaula rule
%   ordrms        - SH information per order (order RMS)
%-----------------------------------------------------------------------
% Complex SH co-efficients to real & vice-versa:
%-----------------------------------------------------------------------
%   cpx2realcov   - Complex to real conversion for covariance matrices of
%                   SH spectrum 
%   cpx2realmat   - Complex to real conversion for a particular SH-degree
%   cpx2realsh    - Complex to real conversion for coefficients
%   LeNorm        - order dependent factors to switch from the geodetic 
%                   normalization of spherical harmonic functions to 
%                   alternative conventions (complex or real representation) 
%   real2cpxcov   - Real to complex conversion for covariance matrices of
%                   SH spectrum 
%   real2cpxmat   - Real to complex conversion for a particular SH-degree
%   real2cpxsh    - Real to complex conversion for coefficients
%-----------------------------------------------------------------------
% Rotation via Wigner-d-function/SO(3)
%-----------------------------------------------------------------------
%  rotate_shc    -  rotation to the spherical harmonic coefficients by 
%                   Wigner-d-function.
%  wigner_all    -  recusive computation of all wigner-d-function for the 
%                   rotation of spherical harmonics
%-----------------------------------------------------------------------
% Slepian functions on the sphere
%-----------------------------------------------------------------------
% evp_gruenbaum  -  solves Gruenbaum's eigenvalue problem for the Slepian 
%                   basis functions which are concentated within a spherical cap
% getSlepianSHC  -  rearranges the coefficients of the Slepian basis functions 
%                    - generated by evp_gruenbaum -  into standard 
%                   |C\S|-format per basis function.
%-----------------------------------------------------------------------
% Normal field:
%-----------------------------------------------------------------------
%   normalklm     - normal field coefficients
%
%-----------------------------------------------------------------------
% Testing:
%-----------------------------------------------------------------------
%   iplmquad      - IPLM by numerical quadrature
%
%-----------------------------------------------------------------------
% Helpers  
%-----------------------------------------------------------------------
%   mex_compile   - contains the command for compiling the mex-functions
%
%-----------------------------------------------------------------------
% Data:
%-----------------------------------------------------------------------
%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% authors: 
%     Nico Sneeuw (NS)
%     Matthias Weigelt (MW)
%     Markus Antoni (MA)
%     Balaji Devaraju (BD)
%     Matthias Roth (MR)
%
% file compilation: 
%   v.2021b: 2021-10-22, MA, update documentation
%                            added: Slepian basis functions within a spherical cap
%   v.2021a: 2021-04-09, MA, remove depricated functions, update documentation
%                            added: Wigner-d-function
%            2015-09-02, MR, add mex_compile
%            2015-07-28, BD, declare some functions deprecated, add replacements
%            2014-10-08, MR, declare some functions deprecated
%   v.2014:  2014-01-15, MR, check all files: revise help texts, beautify code
%   v.2013b: 2013-02-13, MR, consolidation of function names, small bug fixes
%     v.5:   2013-01-18, MA/BD, integration of new tools/update of contents-file
%     v.4:   2012-06-12, MW, added or revision: cs2sc, cs2vec, degcorr, gaussian,
%                            gradpshs, gshsag, gshsptw, nablaPot, normg, pshs, 
%                            sc2cs, tenspshs, updwcont, vec2cs
%     v.3:   2008-08-12, NS, included pointwise shs and numerous other files
%     v.2:   2002-05-22, NS, included ylm, shplot, ylmplot and ancillary files
%     v.1:   2000-10-05, NS, initial version






   
                                       
