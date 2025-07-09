% Tools for processing GRACE spherical harmonic (SH) coefficients
% Requires SHbundle, uberall, FilterBundle and RotationBundle toolboxes
% 
%-----------------------------------------------------------------------
%           Functions for rearranging the SH coefficients
%-----------------------------------------------------------------------
%   klmtsrs2sc  - Rearranges the output of KLMTSRSMAT back in the format provided by READGRC
%   klmtsrsmat  - Rearranges the GRACE SH coefficients in order to perform time-series analysis.
%
%   prpfltr     - Arranges the filter matrices in the required format (required by MASSESTMTR).
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
%-----------------------------------------------------------------------
%                       Time-series analysis
%-----------------------------------------------------------------------
%   lombscargle - Periodogram analysis of time-series (spectrum and space).
%   lssa        - Least-squares spectral analysis of time-series.
%   intannum    - Gives information about longest contiguous part of the time-series.
%   monthfix    - Fills data gaps in the time-series with NaN values.
%
%
%-----------------------------------------------------------------------
%                       Calendar Utilities
%-----------------------------------------------------------------------
%   daysinmonth - Provides the number of days in a month depending on the year.
%   daysinyear  - Provides the number of days in the year.
%   dayofyear   - Provides the day of the year for a given date.
%   doymonth    - Provides the month for a given day of the year
%   doy2cal     - Provides calendar date for a given day of the year
%   grctimetag  - Provides the appropriate time-tag in days for every GRACE monthly solution.
%   isleap      - Determines if a given year is a leap year or not
%   monthnames  - Provides names of months in different formats
%
%
%-----------------------------------------------------------------------
%                   GRACE data specific functions
%-----------------------------------------------------------------------
% I/O functions
%--------------
%   readgrc     - Read GRACE SH coefficients provided in text files by CSR,GFZ and JPL in ICGEM format
%   rdslrc20    - Read SLR C20 values from technical note 7
%   readgrccov  - Read GRACE variance-covariance matrices provided on request by GFZ
%
% Compute from GRACE SH data
%---------------------------
%   dealias     - Tide de-aliasing of GRACE time-series of coefficients as well as grid points
%   detrendgrc  - De-trending the GRACE coefficients as well as grid points
%   getgracemean- Computes the mean of the SH data given a start and end date
%   massestmtr  - Converts GRACE time-series of SH coefficients to spatial grids
%   monthlyres  - Computes the difference between the time-series and its mean-annual cycle.
%   remgracemonths  - Removes GRACE months that are regularized and are flagged '2' by GFZ
%
% Simulating GRACE normal matrices
%---------------------------------
%   NOTE: The following programs require the GRACE L1b position data.
%   blddsgn     - Build the design matrix for computing the GRACE normal matrices 
%   blockinv    - Block inversion of GRACE normal matrices via Schur complement method
%   grcnrml     - Simulates the normal matrix for a GRACE monthly solution.
%   grccvrnc    - Wrapper function for GRCNRML to simulate covariance matrices recursively for the whole time-series.
%   prpcoord    - Prepares GRACE-A GRACE-B positions for GRCNRML
%
%
% Functions for computing catchment specific quantities
%------------------------------------------------------
%   catchagg    - Aggregates mass values over a given catchments.
%   catchleak   - Computes the leakage into all the catchments for a given filter.
%   catchtsrs   - Rearranges the aggregated mass values from MASSESTMTR for performing time-series analysis.
%   cindxcoord  - Provides the coordinates of the pixels of a given catchment.
%   ctchmntcov  - Propagates the GRACE covariance to the catchments.
%   ctchmntdsgn - Provides a spherical harmonic design matrix for the aggregated catchments.
%   fillcatch   - Fill a particular catchment pixels with a given value.
%   findindx    - Retrieves time-series of specified catchments.
%
% Data required for catchment specific quantities
%------------------------------------------------
%   ctchnms.mat - MAT-file containing the catchment names, their area and identification numbers.
%   ctchmntindx.mat     - MAT-file containing the catchment boundaries in a 0.5 x 0.5 degree raster grid.
%   ctchmntindx3.mat    - Same as ctchmntindx.mat but includes Greenland and Antarctica.
%
%-----------------------------------------------------------------------
%		    		Homogeneous isotropic filters
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
%   
%
%-----------------------------------------------------------------------
%                       Visualization tools
%-----------------------------------------------------------------------
%   mapfield    - Maps a grid of mass estimates or other quantities.
%   sctriplot   - Plots the spherical harmonic coefficients as a triangle (SC-format).
%   tsplot      - Plots the time-series from multiple datasets for specified catchments, and also for specified time periods.
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% authors: 
%   Balaji Devaraju (BD)
%   Christof Lorenz (CL)
%   Mohammad J. Tourian (MJT)
%
% file compilation: 
%   v.2013: 2013-02-15, BD, initial set up with the following files
%   v.2017: 2017-05-12, BD, rearranging files for publishing online
