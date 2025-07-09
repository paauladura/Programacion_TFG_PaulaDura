function bf = catchleak(fil,lmax,typ,cmask)

% Computes the leakage of signal into catchments, i.e. leakage-in as is
% commonly referred to in literature
%
% cleak = catchleak(fil,typ,catchmask)
%
% INPUT
% fil 	- 	A [(L+1)^2 * (L+1)^2] filter matrix arranged in order-leading 
%           format, or, a structure with two variables:
%               1. a character array with the filter matrix file;
%               2. path where the file is stored (this order must be
% 				maintained);
%           or, a character array of the filename and the path together, 
%           or, if only the variances are available then provide them
%               in either SC, CS, or order-leading look-up-table formats.
%           Use the first option only when the maximum degree of the
%           spherical harmonic expansion is [lmax <= 70].
% lmax 	-   Maximum degree of the spherical harmonic expansion.
% typ 	-   type of variance-covariance propagation
%               1. 'full' - full spectral variance-covariance matrix is
%                           propagated.
%               2. 'block'- block-diagonal variance-covariance matrix is
%                           propagated.
%               3. 'diag' - only the diagonal of the variance-covariance
%                           matrix is propagated.
% 				4. 'cs'	  - if input is in CS/SC/[l m Clm Slm] formats.
% cmask - 	Catchment mask. It must be a grid with the latitude running
%			from 90 -> -90 and longitude running from 0 -> 360.
%
% OUTPUT
% bf	- 	Leakage of signal into pixels that belong to catchments
% 			stored as a global grid. The pixels outside the catchments
% 			are filled with NaN values (-9999).
%
%------------------------------------------------------------------------
%
% USES 	FilterBundle/
% 		uberall/
% 		SHbundle/gshscovfn
%
%------------------------------------------------------------------------

% Stuttgart, 3 July 2013. Balaji Devaraju

if nargin < 4
	error('All inputs are mandatory')
end

cn = unique(cmask(:));
cn = cn(cn>0);

[r,c] 	 = size(cmask);
dt 		 = pi/r;
[lam,th] = meshgrid(dt/2:dt:(2*pi),(dt/2:dt:pi));

bf = NaN(size(cmask));
ng = struct('size',dt,'length','global');
f  = struct('fil','none');

% Computing cell's Area
theta 	= (0:dt:pi)';
Area	= (cos(theta(1:end-1)) - cos(theta(2:end)))*ones(1,r*2)*(dt^2);

constants

for k = 1:length(cn)
	tmp = find(cmask==cn(k));
	for n = 1:length(tmp)
		b 	= gshscovfn(fil,lmax,[th(tmp(n)) lam(tmp(n))],'none','cell',ng,f,0,typ);
		b 	= b.*Area;
		msk = double(cmask~=cn(k));
		b 	= b.*msk;
		bf(tmp(n)) 	= sum(b(:))/(4*pi);
		% d	= geo2topo([th(:) lam(:)], [th(tmp(n)) lam(tmp(n))]);
		% b	= sptgauss(d*ae,350,1);
		% b 	= b.*Area(:);
		% msk = double(cmask(:)~=cn(k));
		% b 	= b.*msk/(4*pi);
		% bf(tmp(n)) = sum(b);
	end
end
