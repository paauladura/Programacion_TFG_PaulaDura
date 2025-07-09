function h = shplot(field,how,clim,sig,err)

% SHPLOT Spherical Harmonic PLOT produces an SC-triangle plot (logarithmic)
% with colorbar, a degree-RMS plot and an  order-RMS plot. 
%   
% HOW:
%    shplot(field) 
%    h = shplot(field,how,clim,sig,err) 
%
% IN:
%    field ...... may contain signal or error SH coefficients, 
%                 either in SC-triangle format or in CS-square format.
%    how ........ flag, determines how FIELD is plotted. Field itself (0,def.),
%                 as gain wrt. JGM-1s errors (1) or as eigenspectra (2).
%    clim ....... log10 color limits [cmin cmax], optional, may be empty
%    sig ........ flag, determines if the Kaula's rule is plotted.
%    err ........ flag, determines if JGM-1S error curve is plotted.
%
% OUT:
%    h .......... vector of handles of the four axes: triangle-axis,
%                 colorbar, degree-RMS and order-RMS.   
%
% USES: 
%    cs2sc, degreeinfo, shprepare, ordrms, kaula, uberall/trapstrip, 
%    jgm1s.mat
%
% REMARKS: 
%    - how = 1: only useful if FIELD are SD's.  
%    - how = 2: only useful if FIELD are eigenspectra. 
%    - how = 0 or 1: Default sig and err are 1. Otherwise 0. 
%    - how = 2: Condition number plot instead of order-RMS. 

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAGP, TU Munich
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2013-02-13: MR, change function names, brush up comments
%    1998-12-21: NS, updates for Matlab V5: remove scimage, figure position, 
%                    axis
%    1996-09-25: NS, SHPREPARE function used (GAINREF and TRAPSTRIP obsolete)
%    1996-02-13: NS, whole M-file brushed up
%    1994-12-06: NS, how-flag instead of gain-flag (i.e. how=2 option added)
%    1994-10-26: NS, initial version
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

% Set up the figure
clf
set(gcf,'Name','Spherical Harmonic Spectrum',...
        'NumberTitle','off',...
        'Units','centimeters',...
        'Position',[4 5 21 15])
set(gcf,'DefaultAxesUnits','centimeters',...
        'DefaultAxesFontSize',10,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',.5,...
        'DefaultLineLineWidth',3)
%        'DefaultAxesColor',[.8 .8 .8],...
point = get(gcf,'Pointer');
set(gcf,'Pointer','watch')


% Defaults and checking
if nargin < 5, err  = 1;  end		% JGM-1S error curves are plotted
if nargin < 4, sig  = 1;  end		% Kaula's rule is plotted
if nargin < 3, clim = []; end		% can be defined below
if nargin < 2, how  = 0;  end		% FIELD itself is plotted.
if nargin < 1, error('I want input!'), end
if ~any(0:2 == how), error('Invalid HOW flag'), end
if how == 2, sig = 0; err = 0; end	% no Kaula or JGM-1s errors

[rows,cols] = size(field);
lmax        = rows -1;
if rows == cols				% apparently CS-format
    field = cs2sc(field);		% convert into SC-format
elseif cols ~= 2*rows - 1
    error('Input FIELD not in required format!')
end

if err
    errfilename = 'jgm1s.mat'; % Specify path to your JGM-1S file (not provided with SHbundle)
end

% Determine degree and order RMS vectors,
% and their min/max for later plotting purposes.
l = 2:lmax;
m = -lmax:lmax;
drms = degreeinfo(field);			% degree RMS
if how == 2
    orms = max(field)./(10.^min(abs(log10(field))));	% condition numbers
    orms = max(orms,1);
else
    orms = ordrms(field);		% order RMS
end
xmin = min(drms(isfinite(drms))); 	% for determination of 'Xlim' 
xmax = max(drms(isfinite(drms))); 	%    of degree rms plot
ymin = min(orms(isfinite(orms)));		% for determination of 'Ylim' 
ymax = max(orms(isfinite(orms)));		%    of order rms plot


% Now the field is prepared by SHPREPARE for plotting
if how == 0 || how == 1
    field = shprepare(field, how);
else
    field = shprepare(field, 0);	
end

% set empty entries to nan
tr = [flipud(triu(ones(size(field, 1)))), [zeros(1, size(field, 1) - 1); rot90(triu(ones(size(field, 1) - 1)), 2)]];
field(~tr) = nan;

pcolor(m, 0:lmax, field); axis image	% triangle plot
shading flat;
set(gca, 'ydir', 'reverse');            % reverse direction of y-axis
if ~isempty(clim), caxis(clim), end

xlabel('S_{lm} \leftarrow order \rightarrow C_{lm}');
ylabel('degree')
h1 = gca;				% main axes (triangle)

% Now it's time for the colorbar
h2 = colorbar('horiz');			% horizontal color bar with handle
vs = 13 * (lmax+1)/(2*lmax+1);		% vertical size (nearly 6.5)
set(h2,'Position',[2 1 13 .4],'Fontsize',8)
set(h1,'Position',[2 3 13 vs])
if how == 1
    set(get(h2,'Xlabel'),'String','[log10]','Fontsize',9)
elseif how == 0
    set(get(h2,'Xlabel'),'String','Coefficient Size [log10]','Fontsize',9)
elseif how == 2
    set(get(h2,'Xlabel'),'String','Eigenvalue [log10]','Fontsize',9)
end
drawnow, hold on

% Let's do the degree RMS's.
h3 = axes('Position',[16 3 4 vs]);	% handle of degree RMS axes
semilogx(drms,l);                   % degree RMS plot with handle
set(h3,'YDir','reverse')
title('degree RMS')
grid, drawnow, hold on

% Let's do the order RMS's.
h4 = axes('Position',[2 11 13 3]);	% handle of order RMS axes
semilogy(m, orms);      			% order RMS plot with handle
if how == 2
    ylabel('Conditioning')
else
    ylabel('order RMS')
end
grid, drawnow, hold on

% Plot, if requested, additional signal and error curves.
if err 
    load(errfilename);
    djgm = degreeinfo(jgm1s_sd);		% Degree-error RMS JGM-1s
    ojgm = ordrms(jgm1s_sd);		    % Order-error RMS JGM-1
    xmin = min(xmin,min(djgm));
    ymin = min(ymin,min(ojgm));
    hold on
    semilogx(h3, djgm, 2:60, 'r--');   % = axes(h3), semilogx(djgm,2:60,'r--')
    semilogy(h4, -60:60, ojgm, 'r--'); % = axes(h4), semilogy(-60:60,ojgm,'r--'); 
end   
if sig
    dkau = kaula(l);
    okau = ordrms(trapstrip([1e-20 1e-20 dkau]'*ones(1,2*lmax+1)));
    xmin = min(xmin,min(dkau));
    ymin = min(ymin,min(okau));
    hold on
    semilogx(h3, dkau,l,'g'); % axes(h3), semilogx(dkau,l,'g') 
    semilogy(h4, m,okau,'g');    % axes(h4), semilogy(m,okau,'g')
end

% Brush up RMS plots.
xmin = 10^(floor(log10(xmin)));
ymin = 10^(floor(log10(ymin)));
xmax = max(xmax,1e-8);
ymax = max(ymax,1e-8);
set(h3, 'XLim', [xmin xmax], 'YLim', [0 lmax])
set(h4, 'XLim', [-lmax lmax], 'YLim', [ymin ymax])


% finalize plot
axes(h1)				% make the triangle current axes
hold off				% release plot
hh = [h1 h2 h3 h4]'; 
set(hh, 'Units', 'normalized')		% back to normal for scaling purposes
set(gcf, 'Pointer', point)

if nargout == 1, h = hh; end
