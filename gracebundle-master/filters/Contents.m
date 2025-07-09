% Tools for filtering on the sphere , v.2013
% Requires SHbundle, uberall, and RotationBundle toolboxes
% 
%
%-----------------------------------------------------------------------
%		    Homogeneous isotropic filters
%-----------------------------------------------------------------------
%   bttrwrth      - Butterworth filter defined for the spectrum
%   diffusionfil  - Diffusion filter
%   gaussian      - Gaussian filter
%   gaussfltr     - Replacement for gaussian. Uses numerical quadrature
%                   instead of recursion formula
%   pellinen      - Pellinen coefficients of the box-car filter on the sphere
%   spkcosine     - Spectral cosine filter
%   sptbttrwrth   - Spatial Butterworth filter
%   sptcosine     - Spatial cosine filter
%   sptgauss      - Generates the spatial weights of the Gaussian filter directly.
%   vonhann       - von Hann filter (2nd order spatial cosine filter)
%
%-----------------------------------------------------------------------
%                 Inhomogeneous anisotropic filters
%-----------------------------------------------------------------------
%   hannoniso     - Anisotropic Gaussian filter proposed by S-C Han et al 2005
%   fanfilter     - Anisotropic filter named 'Fan' designed by Zhang et al 2009
%   dstrpngmtrx   - Destriping filter of Swenson & Wahr 2006
%
%-----------------------------------------------------------------------
%              Legendre transform: Synthesis & Analysis
%-----------------------------------------------------------------------
%   spt2spk       - Analysis function for the Legendre transform
%   spk2spt       - Synthesis function for the Legendre transform
%
%-----------------------------------------------------------------------
%                  Performance analysis of filters
%-----------------------------------------------------------------------
%   anifperf      - Analyses the performance of inhomgeneous (an)isotropic filters
%   isofperf      - Analyses the performance of homogeneous isotropic filters
%   hisospvar     - Spatial variance of the homogeneous isotropic filters
%   mtfhiso       - Modulation transfer function of the homogeneous isotropic filters
%
%-----------------------------------------------------------------------
%              Format conversion of filter coefficients
%-----------------------------------------------------------------------
%   colombo       - Arranging [l m Clm Slm] in order-leading format
%   colomboQ      - Arranging a filter/covariance matrix in order-leading format
%   degordrngvec  - Arranging [l m Clm Slm] in degree-leading format
%   degordrngQ    - Arranging a filter/covariance matrix in degree-leading format
%   vcm2vec       - Extract the diagonal elements of a degree-leading filter/covariance
%                   matrix and arrange it in [l m Clm Slm] format
%
%-----------------------------------------------------------------------
%                     Rotation of filter matrices
%-----------------------------------------------------------------------
%   rotklm        - Rotating an SH spectrum of a field to a given location
%   rotspkcov     - Rotating the spectrum of a filter to the desired location
%
%-----------------------------------------------------------------------
%                     Power spectrum
%-----------------------------------------------------------------------
%   degvar        - Same as SHbundle/degreeinfo but only provides degree variances
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% authors: 
%     Balaji Devaraju (BD)
%     Nico Sneeuw (NS)
%     Matthias Weigelt (MW)
%
% file compilation: 
%   v.2013: 2013-02-15, BD, initial set up with the following files
%                           anifperf, bttrwrth, colombo, colomboQ, degordrngQ, degordrngvec,
%                           degvar, diffusionfil, dstrpngmtrx, gaussian, hannonsio, isofperf,
%                           pellinen, rotklm, rotspkcov, spk2pt, spkcosine, spt2spk, 
%                           sptbttrwrth, sptcosine, sptgauss, vcm2vec, vonhann
%
%

         
