function Wcs = fanfilter(fil,lmax,mfil)

% FANFILTER computes the spherical harmonic filter coefficients of the Fan
% filter defined in the Zhang et al, 2009. Geophysical Res Lett
%
% fcs = fanfilter
% fcs = fanfilter(fil,lmax,mfil)
%
%
% INPUT
% fil   - Character or structure variable. You can either give the name of a 
%         filter and/or its parameters.
%                   Filter          Parameters 
%               1.  Gauss               cap
%               2.  vonHann             cap
%                   (spatial Cosine filter of order 2)
%               3.  spCos               cap, k
%               4.  Box-car             
%               5.  Butterworth         lc, k
%               6.  Pellinen            cap
%               7.  Diffusion           lc, k
%                   (spectral analog of Gauss filter)
%               8.  Cosine              ls, lc, k
%                   (spectral analog of spatial Cosine filter)
%               9.  spButter            cap, k
%                   (spatial analog of spectral Butterworth filter)
%           cap - filter radius [radians]
%           ls  - start degree [integer]
%           lc  - cut-off degree [integer]
%           k   - order of the filter [integer]
%           For further details on the filter parameters refer to the help text
%           of the individual filter functions.
%           Example: fil = struct('name','Cosine','ls',15,'lc',60,'k',2)
%           Default: fil = struct('name','Gauss','cap',pi/60)
% lmax  - Maximum degree of spherical harmonic expansion [def: 180]
% mfil  - If a different filter for the order direction is desired then provide
%         a filter for it. The format is the same as fil [def: mfil = fil]
%
% OUTPUT
% Wcs   - Filter coefficients given in CS-format
%
%------------------------------------------------------------------------------
% USES
% 	FilterBundle/bttrwrth
% 	             cssc2clm
% 	             degvar
% 	             diffusionfil
% 	             gaussfltr
% 	             pellinen
% 	             spkcosine
% 	             spk2spt
% 	             sptbttrwrth
%                sptcosine
% 	             vonhann
%------------------------------------------------------------------------------

% Initial version: Balaji Devaraju. Stuttgart, 2 April 2014
% Authors: Balaji Devaraju

% Revision History:
%
%------------------------------------------------------------------------------

if nargin == 0
    fil     = 'Gauss';
    lmax    = 180;
    filp    = [];
    mfil    = fil;
    mfilp   = filp;
elseif nargin == 1
    if ischar(fil)
        filp = [];
    elseif isstruct(fil)
        filp    = fil;
        fil     = fil.name;
    else
        error('Variable FIL must be a character string or a structure variable. Please verify.')
    end
    lmax    = 180;
    mfil    = fil;
    mfilp   = filp;
elseif nargin == 2
    if ischar(fil)
        filp = [];
    elseif isstruct(fil)
        filp    = fil;
        fil     = fil.name;
    else
        error('Variable FIL must be a character string or a structure variable. Please verify.')
    end
    mfil    = fil;
    mfilp   = filp;
    if isscalar(lmax)
        lmax = fix(lmax);
    else
        error('LMAX must be a scalar integer.')
    end
elseif nargin == 3
    if ischar(fil)
        filp = [];
    elseif isstruct(fil)
        filp    = fil;
        fil     = fil.name;
    else
        error('Variable FIL must be a character string or a structure variable. Please verify.')
    end
    if ischar(mfil)
        mfilp = [];
    elseif isstruct(mfil)
        mfilp    = mfil;
        mfil     = mfil.name;
    else
        error('Variable FIL must be a character string or a structure variable. Please verify.')
    end
    if isscalar(lmax)
        lmax = fix(lmax);
    else
        error('LMAX must be a scalar integer.')
    end
end

Wl = genfil(fil,lmax,filp);
Wl = Wl/Wl(1);
Wm = genfil(mfil,lmax,mfilp);
Wm = Wm/Wm(1);

tmpl = Wl * ones(1,lmax+1);
tmpm = ones(lmax+1,1) * Wm';
Wcs  = tmpl .* tmpm;
Wcs  = tril(Wcs);
Wcs  = [fliplr(Wcs(:,2:end)) Wcs];
Wcs  = sc2cs(Wcs);

end

function W = genfil(f,L,fp)

% Helper function for generating the desired homogeneous isotropic filter 
% spectrum.


switch f
    case 'Gauss'
        if isempty(fp)
            fp.cap    = pi/36;
        end
        W = gaussfltr(fp.cap,L);
    case 'vonHann'
        if isempty(fp)
           fp.cap = pi/18;
        end
        W = vonhann(fp.cap,L);
    case 'spCos'
        if isempty(fp)
            fp.cap  = pi/18;
            fp.k    = 2;
        end
        W = sptcosine(fp.cap,fp.k,'Gauss',L,true);
    case 'Pellinen'
        if isempty(fp)
            fp.cap    = pi/18;
        end
        W = pellinen(fp.cap,(0:L)');
    case 'Box-car'
        if isempty(fp)
            fp.lc = fix(L/4);
        end
        W = [ones(fp.lc+1,1); zeros(L-fp.lc,1)];
    case 'Butterworth'
        if isempty(fp)
            fc.lc = 15;
            fc.k  = 5;
        end
        W = bttrwrth(fp.lc,fp.k,L);
    case 'Diffusion'
        if isempty(fp)
            fp.lc = 15;
            fp.k  = 1;
        end
        W = diffusionfil(fp.lc,fp.k,L);
    case 'Cosine'
        perf.f = 'Spectral cosine';
        if isempty(fp)
            fp.lc = 15;
            fp.k  = 1;
            if L > 60
                fp.ls = 60;
            else
                fp.ls = L;
            end
        end
        W = spkcosine(fp.lc,fp.k,fp.ls,L);
    case 'spButter'
        if isempty(fp)
            fp.cap  = pi/18;
            fp.k    = 2;
        end
        W = sptbttrwrth(fp.cap,fp.k,'Gauss',L,true);
end

end
