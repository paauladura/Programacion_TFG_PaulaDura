function perf = isofperf(fil,lmax,filp,smpls,fld,nfld,gph)

% ISOFPERF computes certain quantitites that describe the performance of
% the desired homogeneous isotropic filter.
%
% perf = isofperf
% perf = isofperf(fil)
% perf = isofperf(fil,lmax)
% perf = isofperf(fil,lmax,filp)
% perf = isofperf(fil,lmax,filp,smpls)
% perf = isofperf(fil,lmax,filp,smpls,fld)
% perf = isofperf(fil,lmax,filp,smpls,fld,nfld)
% perf = isofperf(fil,lmax,filp,smpls,fld,nfld,gph)
%
% INPUT
% fil   -   A string with the one of the following filter names
%               1. Gauss
%               2. vonHann (spatial Cosine filter of order 2)
%               3. spCosine
%               4. Box-car
%               5. Butterworth
%               6. Pellinen
%               7. Diffusion (spectral analog of Gauss filter)
%               8. Cosine (spectral analog of spatial Cosine filter)
%               9. spButter (spatial analog of spectral Butterworth filter)
%
% lmax  -   Maximum degree of spherical harmonic expansion desired
% filp  -   Structure containing the parameters for the filter provided in
%           'fil'.
%           e.g. Gauss 500 km
%           fil = 'Gauss'
%           filp = struct('cap',500/ae) [ae is the semi-major axis of the
%                   reference ellipsoid].
%
%           e.g. Butterworth
%           fil = 'Butterworth'
%           filp = struct('lc',25,'k',5)
%           'cap' if filter radius in radians
%           'lc' is the cut-off degree and 'k' is the order of the filter.
%           'ls' start degree for cosine taper.
%           Please refer to the help of the different filter functions for 
%           more details.
%
%           e.g. Any filter whose spectrum is known
%           name    = 'other'
%           spk     = Bl (Legendre spectrum of the filter)
%
% smpls -   Sampling in the spatial domain. 
% fld   -   Field for which the quantities have to be computed. This must
%           be given in 'CS', 'SC', or [l m Clm Slm] formats. [optional]
% nfld  -   Noise field in 'CS','SC', or [l m Clm Slm] formats for
%           computing processing gain. [optional]
% gph   -   Toggle switch for plotting the filter curves. [optional] [0 or 1]
%
% OUTPUT
% perf  -   Structure variable containing the performance measures
%           Damping factor, Processing loss, Main-lobe half-width,
%           Main-lobe energy concentration, Processing gain, Spatial
%           leakage, Highest side-lobe level, side-lobe roll-off ratio, 
%           spatial resolution, spectral weights and spatial weights
%--------------------------------------------------------------------------
% See also anifperf, mtfhiso
%--------------------------------------------------------------------------

% USES
% 	FilterBundle/bttrwrth
% 	             degvar
% 	             diffusionfil
% 	             gaussfltr
% 	             pellinen
% 	             spkcosine
% 	             spk2spt
% 	             sptbttrwrth
%                sptcosine
% 	             vonhann
%                hisospvar
%                mtfhiso
%
% 	SHbundle/kaula
%            cssc2clm
%
%--------------------------------------------------------------------------

% Authors:
% Balaji Devaraju. Stuttgart, 23/5/2009

% Revision history
% BD 13/02/2013 - Code brush-up; added comments
% BD 02/04/2014 - Included Spatial Cosine filter and brushed up code.
% BD 08/04/2014 - Included spectral and spatial weights in the output.
% BD 03/04/2014 - Fixed some bugs
%--------------------------------------------------------------------------

if nargin == 0
    lmax = 360;
    filp = [];
    fld = kaula((0:lmax)');
    fld(1:2) = 0;
    fil = 'Gauss';
    smpls = 1000;
    nfld = ones(lmax+1,1);
    gph = false;
elseif nargin == 1
    filp = [];
    lmax = 360;
    fld = kaula((0:lmax)');
    fld(1:2) = 0;
    smpls = 1000;
    nfld = ones(lmax+1,1);
    gph = false;
elseif nargin == 2
    filp = [];
    fld = kaula((0:lmax)');
    fld(1:2) = 0;
    smpls = 1000;
    nfld = ones(lmax+1,1);
    gph = false;
elseif nargin == 3
    fld = kaula((0:lmax)');
    fld(1:2) = 0;
    smpls = 1000;
    nfld = ones(lmax+1,1);
    gph = false;
elseif nargin == 4
    fld = kaula((0:lmax)');
    fld(1:2) = 0;
    nfld = ones(lmax+1,1);
    gph = false;
elseif nargin == 5
    if isempty(smpls)
        smpls = 1000;
    end
    if ~isempty(fld)
        fld = degvar(fld,lmax);
        fld = sqrt(fld(:,2)./(2*fld(:,1) + 1));
    else
        fld = kaula((0:lmax)');
        fld(1:2) = 0;
    end
    nfld = ones(lmax+1,1);
    gph = false;
elseif nargin == 6
    if isempty(smpls)
        smpls = 1000;
    end
    if ~isempty(fld)
        fld = degvar(fld,lmax);
        fld = sqrt(fld(:,2)./(2*fld(:,1) + 1));
    else
        fld = kaula((0:lmax)');
        fld(1:2) = 0;
    end
    nfld = degvar(nfld,lmax);
    nfld = sqrt(nfld(:,2)./(2*nfld(:,1) + 1));
    gph = false;
elseif nargin == 7
    if isempty(smpls)
        smpls = 1000;
    end
    if ~isempty(fld)
        fld = degvar(fld,lmax);
        fld = sqrt(fld(:,2)./(2*fld(:,1) + 1));
    else
        fld = kaula((0:lmax)');
        fld(1:2) = 0;
    end
    nfld = degvar(nfld,lmax);
    nfld = sqrt(nfld(:,2)./(2*nfld(:,1) + 1));
    if gph>=0 && gph<=1
        gph = logical(gph);
    else
        error('GPH variable can take only values 0 and 1. Please check your inputs.')
    end
end


fldnm = fieldnames(filp);

perf.lmax   = lmax;

switch fil
    case 'Gauss'
        perf.fil = 'Gauss';
        if any(strcmp('cap',fldnm))
            perf.cap    = filp.cap;
        else
            perf.cap    = pi/36;
        end
        B = gaussfltr(perf.cap,lmax);
        B = B(1:lmax+1);
    case 'vonHann'
        perf.fil = 'vonHann';
        if any(strcmp('cap',fldnm))
            perf.cap    = filp.cap;
        else
            perf.cap    = pi/18;
        end
        B = vonhann(perf.cap,lmax);
    case 'spCosine'
        perf.fil = 'spCosine';
        if any(strcmp('cap',fldnm))
            perf.cap    = filp.cap;
        else
            perf.cap    = pi/18;
        end
        if any(strcmp('k',fldnm))
            perf.k      = filp.k;
        else
            perf.k      = 2;
        end
        B = sptcosine(perf.cap,perf.k,'Gauss',lmax,true);
    case 'Pellinen'
        perf.fil = 'Pellinen';
        if any(strcmp('cap',fldnm))
            perf.cap    = filp.cap;
        else
            perf.cap    = pi/18;
        end
        B = pellinen(perf.cap,(0:lmax)');
    case 'Box-car'
        perf.fil    = 'Box-car';
        if any(strcmp('lc',fldnm))
            perf.lc = filp.lc;
        else
            perf.lc = 20;
        end
        B = [ones(perf.lc+1,1); zeros(lmax-perf.lc,1)];
    case 'Butterworth'
        perf.fil    = 'Butterworth';
        if any(strcmp('lc',fldnm))
            perf.lc = filp.lc;
        else
            perf.lc = 20;
        end
        if any(strcmp('k',fldnm))
            perf.k      = filp.k;
        else
            perf.k      = 2;
        end
        B = bttrwrth(perf.lc,perf.k,lmax);
    case 'Diffusion'
        perf.fil = 'Diffusion';
        if any(strcmp('lc',fldnm))
            perf.lc = filp.lc;
        else
            perf.lc = 20;
        end
        if any(strcmp('k',fldnm))
            perf.k      = filp.k;
        else
            perf.k      = 1;
        end
        B = diffusionfil(perf.lc,perf.k,lmax);
    case 'Cosine'
        perf.fil = 'Spectral cosine';
        if any(strcmp('lc',fldnm))
            perf.lc = filp.lc;
        else
            if lmax > 60
                perf.lc = 60;
            else
                perf.lc = lmax;
            end
        end
        if any(strcmp('k',fldnm))
            perf.k      = filp.k;
        else
            perf.k      = 1;
        end
        if any(strcmp('k',fldnm))
            perf.ls      = filp.ls;
        else
            perf.ls      = 0;
        end
        B = spkcosine(perf.lc,perf.k,perf.ls,lmax);
    case 'spButter'
        perf.fil = 'Spatial Butterworth';
        if any(strcmp('cap',fldnm))
            perf.cap    = filp.cap;
        else
            perf.cap    = pi/18;
        end
        if any(strcmp('k',fldnm))
            perf.k      = filp.k;
        else
            perf.k      = 2;
        end
        B = sptbttrwrth(perf.cap,perf.k,'Gauss',lmax,true);
    case 'other'
        perf.fil = 'other';
        if ~any(strcmp('spk',fldnm))
            error('Insufficient input arguments')
        else
            B = filp.spk;
        end
    otherwise
        error('String not recognized. Please check the filter name given.')
end

data    = (0:pi/smpls:pi)';
b       = spk2spt(B,lmax,data); % converting spectral windows to spatial kernels

l = (0:lmax)';
m = 2*l + 1; % number of orders in each degree

% Plot graphs of the filters
if gph
    figure
    
    h1 = subplot(2,2,1);
    plot(h1,l,B/B(1),'LineWidth',1)
    axis(h1,[0 200 -0.2 1])
    pbaspect(h1,[1 1 1])
    xlabel(h1,'SH degrees'), ylabel(h1,'Weights')
    grid(h1,'on')
    title(perf.fil)
    
    h2 = subplot(2,2,2); 
    plot(h2,l,20*log10(abs(B/B(1))),'LineWidth',1)
    axis(h2,[0 200 -120 0])
    pbaspect(h2,[1 1 1])
    xlabel(h2,'SH degrees'), ylabel(h2,'Weights [dB]')
    grid(h2,'on')
    
    h3 = subplot(2,2,3); 
    plot(h3,6378.1363*data,b/b(1),'LineWidth',1)
    axis(h3,[0 2000 -0.2 1])
    pbaspect(h3,[1 1 1])
    xlabel(h3,'Spherical distance [km]'), ylabel(h3,'Weights')
    grid(h3,'on')
    
    h4 = subplot(2,2,4); 
    plot(h4,6378.1363*data,20*log10(abs(b/b(1))),'LineWidth',1)
    axis(h4,[0 2000 -120 0])
    pbaspect(h4,[1 1 1])
    xlabel(h4,'Spherical distance [km]'), ylabel(h4,'Weights [dB]')
    grid(h4,'on')
end

perf.M 		= data(b<=0);  % Main-lobe width: First zero-crossing
if isempty(perf.M)
    perf.M = NaN;
else
    perf.M 		= perf.M(1);
end
perf.M 		= perf.M;

if ~isnan(perf.M)
%    perf.sdlb 	= 20*log10(max(abs(b(data>perf.M))/b(1))); % Highest side-lobe
    perf.mlbc 	= sum(b(data<=perf.M).^2 .* sin(data(data<=perf.M)))/sum(b.^2 .* sin(data)); % Mian-lobe concentration
%    perf.rllf 	= 20*log10(abs(b(end)/b(1))) - perf.sdlb; % side-lobe roll-off
else
%    perf.sdlb   = NaN;
    perf.mlbc   = NaN;
%    perf.rllf   = NaN;
end

perf.dmpf 	= sum((B.*fld).^2 .* m)/sum(fld.^2 .* m); % Damping factor
perf.dmpfw 	= sum(B.^2 .* m)/(lmax+1)^2; % Damping factor for field with a white spectrum

perf.ploss 	= 10*log10(1 - perf.dmpf); % Processing loss
perf.plossw = 10*log10(1 - perf.dmpfw); % Processing loss for field with a white spectrum

perf.pgain 	= sum((B.*fld).^2 .* m) * sum(nfld.^2 .* m)/(sum(fld.^2 .* m) * sum((nfld.*B).^2 .* m));
perf.pgain 	= 10*log10(perf.pgain); % processing gain

perf.pgainw = -10*log10(perf.dmpfw); % processing for field with a white spectrum


[perf.spvar,perf.rlen,perf.ormat] = hisospvar(B,lmax);

perf.spvarc 	= sum(b(data<=perf.spvar).^2 .* sin(data(data<=perf.spvar)))/sum(b.^2 .* sin(data)); % Mian-lobe concentration

perf.B      = [l B]; % spectral weights
perf.b      = [data b]; % unnormalized spatial weights


rfil        = struct('name','other','spk',B,'lmax',perf.lmax);
perf.mtf    = mtfhiso(rfil,1.8e4,pi/3600,pi/360,pi/3); % modulation transfer function
if strcmp(perf.fil,'Pellinen') % Special case treatment for Pellinen due to Gibbs effect
    ind             = (perf.mtf(:,1) <= (1.85*perf.cap));
    perf.mtf(ind,2) = 0;
end

ind         = find(perf.mtf(:,2)>1e-8,1,'first');
perf.res    = perf.mtf(ind,1); % resolution of the filter

perf.slkg 	= sum(b(data>perf.res).^2 .* sin(data(data>perf.res)))/sum(b.^2 .* sin(data));% Spatial leakage

if any(strcmp('cap',fldnm))
    rfil = struct('name',perf.fil,'cap',perf.cap);
    if any(strcmp('k',fldnm))
        rfil.k = perf.k;
    end
    perf.mtfnat = mtfhiso(rfil,1.8e4,pi/3600,pi/360,pi/3);

    ind         = find(perf.mtf(:,2)>1e-8,1,'first');
    perf.resnat = perf.mtfnat(ind,1); % resolution of the filter
end

perf.sdlb 	= 20*log10(max(abs(b(data>perf.res))/b(1))); % Highest side-lobe
perf.rllf 	= 20*log10(abs(b(end)/b(1))) - perf.sdlb; % side-lobe roll-off
