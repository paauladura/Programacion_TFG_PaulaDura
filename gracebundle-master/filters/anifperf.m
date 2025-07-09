function [perf,fltr] = anifperf(B,lmax,itype,loc,smpls,fld,nfld,frac,quant)

% ANIFPERF computes the performance of anisotropic filter kernels in the
% spatial domain, i.e., filters that are developed in the spectral domain
% will be analysed for their performance in the spatial domain. This
% function can also be used as for homogeneous isotropic kernels, but will
% be a lot more time consuming. FLTRPERF is specifically written for
% homogeneous isotropic filters and hence it is strongly recommended.
%
% Note: All the calculations are done on the grid points and not in the
%       grid centers. Test calculations have shown hardly any improvement
%       when the calculations were made in grid centers or on a Gaussian
%       grid. Further, a Gaussian grid helps only when calculations are
%       made over the entire range of the latitude. In our case most of the
%       computations are truncated at a certain latitude, and hence, a
%       Gaussian grid does not help the cause.
%
% perf = anifperf
% perf = anifperf(B,lmax,itype)
% perf = anifperf(B,lmax,itype,loc,smpls,fld,nfld,frac,quant)
%
% [perf,fltr] = anifperf(...)
%
% INPUT
% B         -   Spectral matrix of the filter kernel given as a filename,
%               if the fully populated matrix or the order-leading
%               block-diagonal matrix is provided. The matrices must be
%               arranged in order-leading format.
%               If only the diagonal elements of the spectral filter matrix
%               are populated, then the filter can be provided as a CS/SC
%               matrix or as a look-uptable [l m Clm Slm].
% lmax      -   Maximum degree of spherical harmonic expansion. [+ve int]
% itype     -   Input type of the spectral filter matrix
%                   'full'  - Fully populated filter matrix
%                   'block' - Order-leading block-diagonal matrix
%                   'cs'    - CS/SC or [l m Clm Slm] formats
% loc       -   Location of the filter kernel. [co-lat long] [radians]
% smpls     -   Number of samples that have to be used during numerical
%               integration.
% fld       -   For computing damping factor, processing loss, processing
%               gain, processing loss, field-specific main-lobe energy
%               concentration, and field-specific spatial leakage all
%               require the use of a specific field. Must be in CS/SC or as
%               a look-up table [l m Clm Slm].
% nfld      -   Covariance matrix for the provided signal field. It can be
%               a fully populated or order-leading block-diagonal
%               covariance matrix arranged in order-leading format, or if
%               only the diagonal error values are available then they can
%               be provided in CS, SC, or look-up table formats as well.
%               For the fully populated and block-diagonal cases the
%               filename of the covariance matrix along with path of
%               storage can be provided as a character array.The 'nfld'
%               variable is a structure variable with the following values
%                   nfld.cov    - covariance matrix
%                   nfld.itype  - input type 'full', 'block', 'diag'
% frac      -   Fraction for computing main-lobe half-width with different
%               definitions.
% quant     -   Isotropic transfer functions that need to be applied to the
%               fields before the computations
%                   'none' (coefficients define the output), 'geoid',
%                   'potential', 'dg' or 'gravity' (grav. anomaly),
%                   'tr' (grav. disturbance), or 'trr' (2nd rad.
%                   derivative), or 'water' (water equivalent height) or
%                   'smd' (surface mass density).
%                                                   - def: 'none'
%
% OUTPUT
% perf      -   A structure variable containing all the measures.
%                   perf.dmpf   - Damping factor
%                   perf.dmpfw  - Damping factor for a white spectrum
%                   perf.pgain  - Processing gain
%                   perf.ploss  - Processing loss
%                   perf.pwloss - Processing loss for a white spectrum
%                   perf.Mzero  - Main-lobe half-width at zero-crossings
%                   perf.Mpeak  - Main-lobe half-width computed via
%                                   half-width of fraction of peak
%                   perf.Menrgy - Main-lobe half-width computed via
%                                   half-width of fraction of energy
%                   perf.spvar  - Spatial variance
%                   perf.beta   - Main-lobe filter energy concentration
%                   perf.betaf  - Main-lobe signal concentration
%                   perf.slkg   - Spatial filter energy leakage
%                   perf.slkgf  - Spatial signal leakage
%                   perf.hsdlb  - Highest side-lobe level
%                   perf.roll   - Side-lobe roll-off ratio
% fltr      -   The rotated filter along with spherical distance and
%               (180 - azimuth) is output as a structure variable.
%--------------------------------------------------------------------------
% USES gshscovfn clm2sc sc2cs cs2sc cssc2clm eigengrav
%      rotspkcov real2cpxsh rotklm degordrngQ colomboQ gshs
% 	   kaula vcm2vec 
%--------------------------------------------------------------------------
% See also fltrperf
%--------------------------------------------------------------------------

% Balaji Devaraju, 14 January 2011, Stuttgart.
%--------------------------------------------------------------------------

fprintf('--------------------------------------------------------- \n')
fprintf('\tAnisotropic filter performance calculator \n')
fprintf('--------------------------------------------------------- \n\n')
fprintf('Checking input arguments ... ')
if nargin == 0
    lmax     = 50;
    B        = sc2cs(clm2sc(hannoniso(lmax,15,500)));
    itype    = 'cs';
    loc      = [pi, pi]/4;
    smpls    = 360;
    fld      = kaula((0:lmax)');
    fld(1:2) = 0;
    fld      = sc2cs(clm2sc(cssc2clm([(0:lmax)' fld],lmax)));
    nfld     = 1e-16*ones(lmax+1);
    nitype   = 'diag';
    frac     = 0.1;
    quant 	 = ones(lmax+1);
elseif nargin == 3
    if ~ischar(itype)
        error('Input type specification must be a string')
    end
        
    loc      = [pi, pi]/4;
    smpls    = 360;
    fld      = kaula((0:lmax)');
    fld(1:2) = 0;
    fld      = sc2cs(clm2sc(cssc2clm([(0:lmax)' fld],lmax)));
    nfld     = 1e-16*ones(lmax+1);
    nitype   = 'diag';
    frac     = 0.1;
    quant 	 = ones(lmax+1);
elseif nargin == 4
    if ~ischar(itype)
        error('Input type specification must be a string')
    end
    
    if size(loc,2) ~= 2
        error('Location must be a [mx2] array')
    else
        if any(loc(:,1)>pi) && any(loc(:,1)<0)
            error('Co-latitude must be 0<=theta<=180')
        end
    end
    
    smpls    = 360;
    fld      = kaula((0:lmax)');
    fld(1:2) = 0;
    fld      = sc2cs(clm2sc(cssc2clm([(0:lmax)' fld],lmax)));
    nfld     = 1e-16*ones(lmax+1);
    nitype   = 'diag';
    frac     = 0.1;
    quant 	 = ones(lmax+1);
elseif nargin == 5
    if ~ischar(itype)
        error('Input type specification must be a string')
    end
    
    if size(loc,2) ~= 2
        error('Location must be a [mx2] array')
    else
        if any(loc(:,1)>pi) && any(loc(:,1)<0)
            error('Co-latitude must be 0<=theta<=180')
        end
    end
    
    fld      = kaula((0:lmax)');
    fld(1:2) = 0;
    fld      = sc2cs(clm2sc(cssc2clm([(0:lmax)' fld],lmax)));
    nfld     = 1e-16*ones(lmax+1);
    nitype   = 'diag';
    frac     = 0.1;
    quant 	 = ones(lmax+1);
elseif nargin == 6
    if ~ischar(itype)
        error('Input type specification must be a string')
    end
    
    if size(loc,2) ~= 2
        error('Location must be a [mx2] array')
    else
        if any(loc(:,1)>pi) && any(loc(:,1)<0)
            error('Co-latitude must be 0<=theta<=180')
        end
    end
    
    [r,c] = size(fld);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        fld = sc2cs(clm2sc(fld));
    elseif isequal(c,(2*lmax+1))
        fld = sc2cs(fld);
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    
    nfld    = 1e-16*ones(lmax+1);
    frac    = 0.1;
    quant   = ones(lmax+1);
elseif nargin == 7
    if ~ischar(itype)
        error('Input type specification must be a string')
    end
    
    if size(loc,2) ~= 2
        error('Location must be a [mx2] array')
    else
        if any(loc(:,1)>pi) && any(loc(:,1)<0)
            error('Co-latitude must be 0<=theta<=180')
        end
    end
    
    [r,c] = size(fld);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        fld = sc2cs(clm2sc(fld));
    elseif isequal(c,(2*lmax+1))
        fld = sc2cs(fld);
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    
    tmp     = fieldnames(nfld);
    nitype  = nfld.(tmp{2});
    tmp     = nfld.(tmp{1});
    if ischar(tmp)
        nfld = load(tmp);
        tmp  = fieldnames(nfld);
        nfld = nfld.(tmp{1});
    else
        nfld = tmp;
    end
    
    [r,c] = size(nfld);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        nfld = sc2cs(clm2sc(nfld));
    elseif isequal(c,(2*lmax+1))
        nfld = sc2cs(nfld);
    elseif isequal(r,c) && isequal(r, (lmax+1)^2)
        if strcmp(nitype,'diag')
            nfld = degordrngQ(nfld,lmax);
            nfld = sc2cs(clm2sc(vcm2vec(nfld,lmax)));
        end
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    frac  = 0.1;
    quant = ones(lmax+1);
elseif nargin == 8
    if ~ischar(itype)
        error('Input type specification must be a string')
    end
    
    if size(loc,2) ~= 2
        error('Location must be a [mx2] array')
    else
        if any(loc(:,1)>pi) && any(loc(:,1)<0)
            error('Co-latitude must be 0<=theta<=180')
        end
    end
    
    [r,c] = size(fld);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        fld = sc2cs(clm2sc(fld));
    elseif isequal(c,(2*lmax+1))
        fld = sc2cs(fld);
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    
    tmp     = fieldnames(nfld);
    nitype  = nfld.(tmp{2});
    tmp     = nfld.(tmp{1});
    if ischar(tmp)
        nfld = load(tmp);
        tmp  = fieldnames(nfld);
        nfld = nfld.(tmp{1});
    else
        nfld = tmp;
    end
    
    [r,c] = size(nfld);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        nfld = sc2cs(clm2sc(nfld));
    elseif isequal(c,(2*lmax+1))
        nfld = sc2cs(nfld);
    elseif isequal(r,c) && isequal(r, (lmax+1)^2)
        if strcmp(nitype,'diag')
            nfld = degordrngQ(nfld,lmax);
            nfld = sc2cs(clm2sc(vcm2vec(nfld,lmax)));
        end
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    
    if frac>1 || frac<=0
        error('Fraction must be less than 1 and greater than zero.')
    end
    quant = ones(lmax+1);
elseif nargin == 9
    if ~ischar(itype)
        error('Input type specification must be a string')
    end
    
    if size(loc,2) ~= 2
        error('Location must be a [mx2] array')
    else
        if any(loc(:,1)>pi) && any(loc(:,1)<0)
            error('Co-latitude must be 0<=theta<=180')
        end
    end
    
    [r,c] = size(fld);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        fld = sc2cs(clm2sc(fld));
    elseif isequal(c,(2*lmax+1))
        fld = sc2cs(fld);
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    
    tmp     = fieldnames(nfld);
    nitype  = nfld.(tmp{2});
    tmp     = nfld.(tmp{1});
    if ischar(tmp)
        nfld = load(tmp);
        tmp  = fieldnames(nfld);
        nfld = nfld.(tmp{1});
    else
        nfld = tmp;
    end
    
    [r,c] = size(nfld);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        nfld = sc2cs(clm2sc(nfld));
    elseif isequal(c,(2*lmax+1))
        nfld = sc2cs(nfld);
    elseif isequal(r,c) && isequal(r, (lmax+1)^2)
        if strcmp(nitype,'diag')
            nfld = degordrngQ(nfld,lmax);
            nfld = sc2cs(clm2sc(vcm2vec(nfld,lmax)));
        end
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    
    if frac>1 || frac<=0
        error('Fraction must be less than 1 and greater than zero.')
    end
    quant = sc2cs(clm2sc(cssc2clm([(0:lmax)' eigengrav((0:lmax)',quant,0)],lmax)));
elseif nargin > 0 && nargin < 3
    error('Insufficient input arguments')
end

fprintf('complete \n\n')

% Indices of starting and ending elements of cosine and sine components
cb = 1;
ce = sum(1:lmax+1);
sb = ce+1;
se = (lmax+1)^2;

% nsph = 180; % Number of latitude samples for computing the performance measures

perf.dmpf         = 0;
perf.dmpfw        = 0;
perf.pgain        = 0;
perf.ploss        = 0;
perf.plossw       = 0;

perf.Mzero        = cell(size(loc,1),1);
perf.Mpeak        = cell(size(loc,1),1);
perf.Menrgy       = cell(size(loc,1),1);
perf.spvar        = cell(size(loc,1),1);
perf.spvarell     = cell(size(loc,1),1);
perf.rlen         = cell(size(loc,1),1);
perf.hsdlb        = cell(size(loc,1),1);
perf.roll         = cell(size(loc,1),1);

perf.beta.Mzero   = cell(size(loc,1),1);
perf.beta.Mpeak   = cell(size(loc,1),1);
perf.beta.spvar   = cell(size(loc,1),1);

perf.cumenergy    = cell(size(loc,1),1);
%perf.betaf.Mzero  = cell(size(loc,1),1);
%perf.betaf.Mpeak  = cell(size(loc,1),1);
%perf.betaf.Menrgy = cell(size(loc,1),1);
%perf.betaf.spvar  = cell(size(loc,1),1);
%
%perf.slkg.Mzero   = cell(size(loc,1),1);
%perf.slkg.Mpeak   = cell(size(loc,1),1);
%perf.slkg.spvar   = cell(size(loc,1),1);
%
%perf.slkgf.Mzero  = cell(size(loc,1),1);
%perf.slkgf.Mpeak  = cell(size(loc,1),1);
%perf.slkgf.Menrgy = cell(size(loc,1),1);
%perf.slkgf.spvar  = cell(size(loc,1),1);
%perf.slkgf.gam    = cell(size(loc,1),1);
%perf.slkgf.rlen   = cell(size(loc,1),1);
%perf.slkgf.ormat  = cell(size(loc,1),1);

fltr.b            = cell(size(loc,1),1);
fltr.clat         = zeros(smpls+1,1);
fltr.long         = zeros(2*smpls,1);

%--------------------------------------------------------------------------
% Calculations for filters of the form B_lm_nk provided in order-leading
% format with (lmax+1)^4 elements
%--------------------------------------------------------------------------
% Computing the design matrices for all the locations
[ylmc,ylms] = ylmdsgn(loc,lmax);

if strcmp(itype,'full') || strcmp(itype,'block')
    if ischar(B)
        B   = load(B);
        tmp = fieldnames(B);
        B   = B.(tmp{1});
    end
    [r,c] = size(B);
    if r~=c || r~=(lmax+1)^2
        error('Filter matrix (B) is not in the required dimensions')
    end
    
    %--------------------------------------------------------------------------
    % Calculating global performance measures
    %--------------------------------------------------------------------------
    
    fprintf('-------------------------------------------------------\n')
    fprintf('\tCalculating global performance measures \n')
    fprintf('-------------------------------------------------------\n')

    fld = fld.*quant;
    fprintf('Damping factor ... ')
    tmp       = cssc2clm(fld,lmax);
    tmp       = [tmp(:,3); tmp(lmax+2:end,4)];
    tmp       = B*tmp;
    perf.dmpf = sum(tmp.^2)/sum(sum(fld.^2));
    fprintf('complete \n')
    
    fprintf('Dampfing factor for a white spectrum ... ')
    perf.dmpfw = sum((sum(B,2)).^2)/(lmax+1)^2;
    fprintf('complete \n')
    
    fprintf('Processing gain ... ')
    if strcmp(nitype,'full') || strcmp(nitype,'block')
        qtmp = cssc2clm(quant,lmax);
        qtmp = [qtmp(:,3); qtmp(lmax+2:end,4)];
        nfld = (qtmp*qtmp').*nfld;
        ntmp = B*nfld;
        ntmp = ntmp*B';
        ntmp = sum(diag(nfld))/sum(diag(ntmp));
    elseif strcmp(nitype,'diag')
        ntmp = B.*B;
        nfld = cssc2clm(nfld,lmax);
        nfld = [nfld(:,3); nfld(lmax+2:end,4)];
        ntmp = ntmp*nfld;
        ntmp = sum(nfld)/sum(ntmp);
    end
    perf.pgain = perf.dmpf*ntmp;
    fprintf('complete \n')
    
    fprintf('Processing loss ... ')
    perf.ploss = 1 - perf.dmpf;
    fprintf('complete \n')
    
    fprintf('Processing loss for a white spectrum ... ')
    perf.plossw = 1- perf.dmpfw;
    fprintf('complete \n\n')
    
    %----------------------------------------------------------------------
    % Calculating local performance measures
    %----------------------------------------------------------------------

    for m = 1:size(loc,1)
        strng = ['Co-lat. ', num2str(loc(m,1)*180/pi),', Long. ', num2str(loc(m,2)*180/pi)];
        fprintf('--------------------------------------------------------------------------\n')
        fprintf('\tCalculating local performance measures for %s \n', strng)
        fprintf('--------------------------------------------------------------------------\n\n')
        
        fprintf('Computing Y*B ...')

        tmp                 = cssc2clm(zeros(lmax+1),lmax);
        if strcmp(itype,'full')
            tmp(:,3)            = ylmc(m,:)*B(cb:ce,cb:ce) + ylms(m,:)*B(sb:se,cb:ce);
            tmp(lmax+2:end,4)   = ylmc(m,:)*B(cb:ce,sb:se) + ylms(m,:)*B(sb:se,sb:se);
        elseif strcmp(itype,'block')
            tmp(:,3)            = ylmc(m,:)*B(cb:ce,cb:ce);
            tmp(lmax+2:end,4)   = ylms(m,:)*B(sb:se,sb:se);
        end

        Bsc = clm2sc(tmp);
        fprintf('complete\n')

        fprintf('Rotating the filter ... ')
        %Bl      = degordrngQ(B,lmax);
        %C       = cpx2realcov(lmax);
        %Bl      = C'*Bl;
        %Bl      = Bl*C;
        %Bl      = rotspkcov(Bl,lmax,loc(m,1),loc(m,2),'full');
        %Bl      = C*Bl;
        %Bl      = Bl*C';
        %Bl      = real(Bl);
        %Bl      = colomboQ(Bl,lmax);
        
        for k = 1:lmax
            C   = cpx2realmat(k);
            D   = diag(exp(1i*(-k:k)'*loc(m,2)))* dlmk(k,-loc(m,1)); % Computation of rotation matrix
            D   = real(C*D*C');
            idx = lmax+1-k:lmax+1+k;
            Bsc(k+1,idx) = Bsc(k+1,idx)*D;
        end

        fprintf('complete \n')
        
        fprintf('Propagating the filter to the spatial domain ... ')
        
        [b,t,l]     = gshs(Bsc,'none','mesh',smpls,0,0);
        fltr.b{m}   = b;
        
        fprintf('done\n')

        % [b,t,l] = gshscovfn(Bl,lmax,loc(m,:),'none','mesh',struct('size',pi/smpls,'length','global'),[],0,'full');
        % if l(end) == 2*pi,
        %     b = b(:,1:end-1);
        %     l = l(1:end-1);
        % end
        % Bl = [];

        fprintf('\n')
        
        fprintf('Rotating the field ... ')
        klm = real2cpxsh(fld,lmax);
        klm = rotklm(klm,lmax,loc(m,1),loc(m,2));
		klm	= cpx2realsh(klm,lmax);
        rf  = gshs(klm,'none','mesh',smpls,0,0);
        fprintf('complete \n\n')
        
        fprintf('Main-lobe half-width ... ')
        perf.Mzero{m} = [l(:) ones(size(l(:)))*NaN];
        idx           = (b<=0);
        for k = 1:5, 
            idx = cumsum(idx); 
        end
        [r,cols]               = find(idx==1);
        perf.Mzero{m}(cols,2)  = t(r);
        fprintf('zero-crossings \n')
        %
        fprintf('Main-lobe half-width ... ')
        idx = (b<=frac*max(b(:)));
        for k = 1:5, 
            idx = cumsum(idx); 
        end
        [r,cols]      = find(idx==1);
        perf.Mpeak{m} = [l(cols)' t(r)];
        fprintf('fraction of the peak \n')
        %
        fprintf('Main-lobe half-width ... ')
        dpsi            = (pi/smpls)*sin(t)*ones(1,length(l));
        Eb              = cumsum(b.^2 .* dpsi);
        totE            = sum(Eb(end,:))*pi/smpls;
        Sb              = cumsum(b.*rf.*dpsi);
        totS            = sum(Sb(end,:))*pi/smpls;
        idx             = bsxfun(@ge, Eb, Eb(end,:)*(1-frac));
        idx             = cumsum(idx);
        [r,cols]        = find(idx==1);
        perf.Menrgy{m}  = [l(cols)' t(r)];
        
        for k = 1:size(Eb,2)
            idx                 = (Eb(:,k) >= Eb(end,k)*(1-frac));
            idx                 = find(idx,1,'first');
            perf.Menrgy{m}(k,2) = t(idx);
        end
        fprintf('fraction of the total energy \n')
        
        fprintf('Spatial variance of the smoothing kernel ... ')
        bsqr    = b.^2 .* dpsi/totE;
        
        % Calculation of mean of the rotated cartesian coordinates
        ups     = sin(t)*cos(l);
        nu      = sin(t)*sin(l);
        zeta    = cos(t)*ones(size(l));

        upsbar  = sum(ups(:).*bsqr(:)*pi/smpls);
        nubar   = sum(nu(:).*bsqr(:)*pi/smpls);
        zetabar = sum(zeta(:).*bsqr(:)*pi/smpls);

        % Calculation of resultant length
        perf.rlen{m}  = sqrt(upsbar^2 + nubar^2 + zetabar^2);

        %Calculation of the variances of the rotated cartesian coordinates
        sigups  = sum(ups(:).^2  .* bsqr(:) * pi/smpls);
        signu   = sum(nu(:).^2   .* bsqr(:) * pi/smpls);
        sigzeta = sum(zeta(:).^2 .* bsqr(:) * pi/smpls);

        sigupsnu    = sum(ups(:) .* nu(:)   .* bsqr(:) * pi/smpls);
        sigupszeta  = sum(ups(:) .* zeta(:) .* bsqr(:) * pi/smpls);
        signuzeta   = sum(nu(:)  .* zeta(:) .* bsqr(:) * pi/smpls);

        % Calculation of orientation matrix
        T   = [sigups sigupsnu sigupszeta; ...
                sigupsnu signu signuzeta; ...
                sigupszeta signuzeta sigzeta];

        % Calculation of spatial variance
        gam = 0.5 * atan2(2*sigupsnu,(sigups - signu));

        R   = [cos(gam) -sin(gam) 0; ...
                sin(gam) cos(gam) 0; ...
                0 0 1];

        TT  = R' * T * R;

        perf.gam{m}     = gam;
        perf.ormat{m}   = T;
        perf.spvar{m}   = asin(sqrt([TT(1,1) TT(2,2)]));
        perf.spvarell{m}= [l(:) asin(sqrt(TT(1,1)*TT(2,2))./sqrt(TT(2,2)*cos(l(:)).^2 + TT(1,1)*sin(l(:)).^2))];

        fprintf('complete \n')
        
        fprintf('Main-lobe concentration ... ')
        perf.beta.Mzero{m}   = perf.Mzero{m};
        perf.beta.Mpeak{m}   = perf.Mpeak{m};
        perf.beta.spvar{m}   = perf.Mpeak{m};
        
        %perf.betaf.Mzero{m}  = perf.Mzero{m};
        %perf.betaf.Mpeak{m}  = perf.Mpeak{m};
        %perf.betaf.Menrgy{m} = perf.Menrgy{m};
        %perf.betaf.spvar{m}  = perf.Mpeak{m};
        %
        %perf.slkg.Mzero{m}   = perf.Mzero{m};
        %perf.slkg.Mpeak{m}   = perf.Mpeak{m};
        %perf.slkg.spvar{m}   = perf.Mpeak{m};
        %
        %perf.slkgf.Mzero{m}  = perf.Mzero{m};
        %perf.slkgf.Mpeak{m}  = perf.Mpeak{m};
        %perf.slkgf.Menrgy{m} = perf.Menrgy{m};
        %perf.slkgf.spvar{m}  = perf.Mpeak{m};
        
        for k = 1:size(Eb,2)
            % Filter energy
            perf.beta.Mzero{m}(k,2)   = Eb(perf.Mzero{m}(k,2) == t,k)/totE;
            %perf.slkg.Mzero{m}(k,2)   = Eb(end,k)/totE - perf.beta.Mzero{m}(k,2);
            perf.beta.Mpeak{m}(k,2)   = Eb(perf.Mpeak{m}(k,2) == t,k)/totE;
            %perf.slkg.Mpeak{m}(k,2)   = Eb(end,k)/totE - perf.beta.Mpeak{m}(k,2);
            perf.beta.spvar{m}(k,2)   = Eb(find(t <= perf.spvarell{m}(k,2),1,'last'),k)/totE;
            %perf.slkg.spvar{m}(k,2)   = Eb(end,k)/totE - perf.beta.spvar{m}(k,2);
            
            % Signal
            %perf.betaf.Mzero{m}(k,2)  = Sb(perf.Mzero{m}(k,2) == t,k)/totS;
            %perf.slkgf.Mzero{m}(k,2)  = Sb(end,k)/totS - perf.betaf.Mzero{m}(k,2);
            %perf.betaf.Mpeak{m}(k,2)  = Sb(perf.Mpeak{m}(k,2) == t,k)/totS;
            %perf.slkgf.Mpeak{m}(k,2)  = Sb(end,k)/totS - perf.betaf.Mpeak{m}(k,2);
            %perf.betaf.Menrgy{m}(k,2) = Sb(perf.Menrgy{m}(k,2) == t,k)/totS;
            %perf.slkgf.Menrgy{m}(k,2) = Sb(end,k)/totS - perf.betaf.Menrgy{m}(k,2);
            %perf.betaf.spvar{m}(k,2)  = Sb(find(t <= perf.spvar{m,1},1,'last'),k)/totS;
            %perf.slkgf.spvar{m}(k,2)  = Sb(end,k)/totS - perf.betaf.spvar{m}(k,2);
        end
        fprintf('complete \n')
        
        %fprintf('Spatial leakage ... ')
        %
        %fprintf('complete \n')
        
        fprintf('Cumulative energy ... ')
        tmp = sum((b.^2)*pi/smpls,2);
        tmp = flipud(cumsum(flipud(sin(t).*tmp*pi/smpls)));
        tmp = bsxfun(@plus,b.^2,[tmp(2:end);0]);
        
        perf.cumenergy{m} = tmp/totE;
        fprintf('complete \n')
        
        fprintf('Highest side-lobe level ... ')
        perf.hsdlb{m}      = [l(:) l(:)];
        perf.hsdlb{m}(:,2) = 20*log10(abs(min(b)./max(b)));
        fprintf('complete \n')
        
        fprintf('Side-lobe roll-off ratio ... ')
        perf.roll{m}   = [l(:) l(:)];
        perf.roll{m}(:,2) = 20*log10(abs(b(end,:)./max(b)));
        perf.roll{m}(:,2) = perf.roll{m}(:,2) - perf.hsdlb{m}(:,2);
        fprintf('complete \n')

    end
    fltr.clat = t;
    fltr.long = l;
    
    
%----------------------------------------------------------------------
% Calculations for filters of the form B_lm_lm provided in CS-format
%----------------------------------------------------------------------
elseif strcmp(itype,'cs')
    [r,c] = size(B);
    if isequal(c,4) && isequal(r, sum(1:lmax+1))
        B = sc2cs(clm2sc(B));
    elseif isequal(c,(2*lmax+1))
        B = sc2cs(B);
    elseif ~isequal(r,c)
        error('Filter matrix (B) not in the specified format')
    end
    
    %----------------------------------------------------------------------
    % Calculating global performance measures
    %----------------------------------------------------------------------
    
    fld = fld.*quant;
    fprintf('Damping factor ... ')
    tmp       = B.*fld;
    perf.dmpf = sum(tmp(:).^2)/sum(fld(:).^2);
    fprintf('complete \n')
    
    fprintf('Dampfing factor for a white spectrum ... ')
    perf.dmpfw = sum(B(:).^2)/(lmax+1)^2;
    fprintf('complete \n')
    
    fprintf('Processing gain ... ')
    if strcmp(nitype,'full') || strcmp(nitype,'block')
        qtmp = cssc2clm(B.*quant,lmax);
        qtmp = [qtmp(:,3); qtmp(lmax+2:end,4)];
        ntmp = (qtmp*qtmp').*nfld;
        qtmp = cssc2clm(quant,lmax);
        qtmp = [qtmp(:,3); qtmp(lmax+2:end,4)];
        nfld = (qtmp*qtmp').*nfld;
        ntmp = sum(diag(nfld))/sum(diag(ntmp));
    elseif strcmp(nitype,'diag')
        ntmp = (B.*quant).^2 .* nfld;
        nfld = nfld.* quant.^2;
        ntmp = sum(nfld(:))/sum(ntmp(:));
    end
    perf.pgain = perf.dmpf*ntmp;
    fprintf('complete \n')
    
    fprintf('Processing loss ... ')
    perf.ploss = 1 - perf.dmpf;
    fprintf('complete \n')
    
    fprintf('Processing loss for a white spectrum ... ')
    perf.plossw = 1- perf.dmpfw;
    fprintf('complete \n\n')
    
    %----------------------------------------------------------------------
    % Calculating local performance measures
    %----------------------------------------------------------------------
    
    for m = 1:size(loc,1)
        strng = ['Co-lat. ', num2str(loc(m,1)),', Long. ', num2str(loc(m,2))];
        fprintf('------------------------------------------------------------------\n')
        fprintf('Calculating local performance measures for %s \n', strng)
        fprintf('------------------------------------------------------------------\n\n')
        
        fprintf('Rotating the filter ... ')
        C       = cpx2realcov(lmax);
        Bl      = rotspkcov(B,lmax,loc(m,1),loc(m,2),'cs');
        Bl      = real(C*Bl*C');
        Bl      = colomboQ(Bl,lmax);
        fprintf('complete \n')
        
        fprintf('Propagating the filter to the spatial domain ... \n')
        tmp                 = cssc2clm(zeros(lmax+1),lmax);
        tmp(:,3)            = ylmc(m,:)*Bl(cb:ce,cb:ce) + ylms(m,:)*Bl(sb:se,cb:ce);
        tmp(lmax+2:end,4)   = ylmc(m,:)*Bl(cb:ce,sb:se) + ylms(m,:)*Bl(sb:se,sb:se);

        Bl = clm2sc(tmp);

        [b,t,l] = gshs(Bl,'none','mesh',smpls,0,0);

        % [b,t,l] = gshscovfn(Bl,lmax,loc(m,:),'none','mesh',struct('size',pi/smpls,'length','global'),[],0,'full');
        % if l(end) == 2*pi,
        %     b       = b(:,1:end-1);
        %     l       = l(1:end-1);
        % end
        % Bl      = [];
        
        fprintf('\n')
        
        fprintf('Rotating the field ...')
        klm     = real2cpxsh(fld,lmax);
        klm     = rotklm(klm,lmax,loc(m,1),loc(m,2));
		klm		= cpx2realsh(klm,lmax);
        rf      = gshs(klm,'none','mesh',smpls,0,0);
        fprintf('complete \n\n')
        
        fprintf('Main-lobe half-width ... ')
        perf.Mzero{m} = [l(:) ones(size(l(:)))*180];
        idx           = (b<=0);
        for k = 1:5, 
            idx = cumsum(idx); 
        end
        [r,cols]               = find(idx==1);
        perf.Mzero{m}(cols,2)  = t(r);
        fprintf('zero-crossings \n')
        %
        fprintf('Main-lobe half-width ... ')
        idx           = (b<=frac*max(b(:)));
        for k = 1:5, 
            idx = cumsum(idx); 
        end
        [r,cols]      = find(idx==1);
        perf.Mpeak{m} = [l(cols)' t(r)];
        fprintf('fraction of the peak \n')
        %
        fprintf('Main-lobe half-width ... ')
        dpsi           = repmat(sin(t)*(pi/smpls),1,length(l));
        Eb             = cumsum(b.^2 .* dpsi);
        totE           = sum(Eb(end,:))*pi/smpls;
        Sb             = cumsum(b.*rf.*dpsi);
        totS           = sum(Sb(end,:))*pi/smpls;
        perf.Menrgy{m} = [l(:) l(:)];
        for k = 1:size(Eb,2)
            idx                 = (Eb(:,k) >= Eb(end,k)*(1-frac));
            idx                 = find(idx,1,'first');
            perf.Menrgy{m}(k,2) = t(idx);
        end
        fprintf('fraction of the total energy \n')
        
        fprintf('Spatial variance of the smoothing kernel ... ')
        bsqr            = b.^2 .* (dpsi .* repmat(t.^2,1,length(l))) / totE;
        
        % Calculation of mean of the rotated cartesian coordinates
        ups     = sin(t)*cos(l);
        nu      = sin(t)*sin(l);
        zeta    = cos(t)*ones(size(l));

        upsbar  = sum(ups(:).*bsqr(:)*pi/smpls);
        nubar   = sum(nu(:).*bsqr(:)*pi/smpls);
        zetabar = sum(zeta(:).*bsqr(:)*pi/smpls);

        % Calculation of resultant length
        perf.rlen{m}  = sqrt(upsbar^2 + nubar^2 + zetabar^2);

        %Calculation of the variances of the rotated cartesian coordinates
        sigups  = sum(ups(:).^2  .* bsqr(:) * pi/smpls);
        signu   = sum(nu(:).^2   .* bsqr(:) * pi/smpls);
        sigzeta = sum(zeta(:).^2 .* bsqr(:) * pi/smpls);

        sigupsnu    = sum(ups(:) .* nu(:)   .* bsqr(:) * pi/smpls);
        sigupszeta  = sum(ups(:) .* zeta(:) .* bsqr(:) * pi/smpls);
        signuzeta   = sum(nu(:)  .* zeta(:) .* bsqr(:) * pi/smpls);

        % Calculation of orientation matrix
        T   = [sigups sigupsnu sigupszeta; ...
                sigupsnu signu signuzeta; ...
                sigupszeta signuzeta sigzeta];

        % Calculation of spatial variance
        gam = 0.5 * atan2(2*sigupsnu,(sigups - signu));

        R   = [cos(gam) -sin(gam) 0; ...
                sin(gam) cos(gam) 0; ...
                0 0 1];

        TT  = R' * T * R;

        perf.gam{m}     = gam;
        perf.ormat{m}   = T;
        perf.spvar{m}   = asin([TT(1,1) TT(2,2)]);

        % perf.spvar{m,2} = [l(:) sum(bsqr)'];
        % perf.spvar{m,1} = sqrt(sum(bsqr(:))*pi/smpls);

        fprintf('complete \n')
        
        %fprintf('Main-lobe concentration ... ')
        %perf.beta.Mzero{m}   = perf.Mzero{m};
        %perf.beta.Mpeak{m}   = perf.Mpeak{m};
        %perf.beta.spvar{m}   = perf.Mpeak{m};
        %
        %perf.betaf.Mzero{m}  = perf.Mzero{m};
        %perf.betaf.Mpeak{m}  = perf.Mpeak{m};
        %perf.betaf.Menrgy{m} = perf.Menrgy{m};
        %perf.betaf.spvar{m}  = perf.Mpeak{m};
        %
        %perf.slkg.Mzero{m}   = perf.Mzero{m};
        %perf.slkg.Mpeak{m}   = perf.Mpeak{m};
        %perf.slkg.spvar{m}   = perf.Mpeak{m};
        %
        %perf.slkgf.Mzero{m}  = perf.Mzero{m};
        %perf.slkgf.Mpeak{m}  = perf.Mpeak{m};
        %perf.slkgf.Menrgy{m} = perf.Menrgy{m};
        %perf.slkgf.spvar{m}  = perf.Mpeak{m};
        %
        %for k = 1:size(Eb,2)
        %    % Filter energy
        %    perf.beta.Mzero{m}(k,2)   = Eb(perf.Mzero{m}(k,2) == t,k)/totE;
        %    perf.slkg.Mzero{m}(k,2)   = Eb(end,k)/totE - perf.beta.Mzero{m}(k,2);
        %    perf.beta.Mpeak{m}(k,2)   = Eb(perf.Mpeak{m}(k,2) == t,k)/totE;
        %    perf.slkg.Mpeak{m}(k,2)   = Eb(end,k)/totE - perf.beta.Mpeak{m}(k,2);
        %    perf.beta.spvar{m}(k,2)   = Eb(find(t <= perf.spvar{m,1},1,'last'),k)/totE;
        %    perf.slkg.spvar{m}(k,2)   = Eb(end,k)/totE - perf.beta.spvar{m}(k,2);
        %    
        %    % Signal
        %    perf.betaf.Mzero{m}(k,2)  = Sb(perf.Mzero{m}(k,2) == t,k)/totS;
        %    perf.slkgf.Mzero{m}(k,2)  = Sb(end,k)/totS - perf.betaf.Mzero{m}(k,2);
        %    perf.betaf.Mpeak{m}(k,2)  = Sb(perf.Mpeak{m}(k,2) == t,k)/totS;
        %    perf.slkgf.Mpeak{m}(k,2)  = Sb(end,k)/totS - perf.betaf.Mpeak{m}(k,2);
        %    perf.betaf.Menrgy{m}(k,2) = Sb(perf.Menrgy{m}(k,2) == t,k)/totS;
        %    perf.slkgf.Menrgy{m}(k,2) = Sb(end,k)/totS - perf.betaf.Menrgy{m}(k,2);
        %    perf.betaf.spvar{m}(k,2)  = Sb(find(t <= perf.spvar{m,1},1,'last'),k)/totS;
        %    perf.slkgf.spvar{m}(k,2)  = Sb(end,k)/totS - perf.betaf.spvar{m}(k,2);
        %end
        %fprintf('complete \n')
        
        fprintf('Spatial leakage ... ')
        
        fprintf('complete \n')
        
        fprintf('Highest side-lobe level ... ')
        perf.hsdlb{m}      = [l(:) l(:)];
        perf.hsdlb{m}(:,2) = 20*log10(abs(min(b)./max(b)));
        fprintf('complete \n')
        
        fprintf('Side-lobe roll-off ratio ... ')
        perf.roll{m}   = [l(:) l(:)];
        perf.roll{m}(:,2) = 20*log10(abs(b(end,:)./max(b)));
        perf.roll{m}(:,2) = perf.roll{m}(:,2) - perf.hsdlb{m}(:,2);
        fprintf('complete \n')
        
        fltr.b{m} = b;
    end
    fltr.clat = t;
    fltr.long = l;
else
    error('Unknown input type for the filter matrix')
end
