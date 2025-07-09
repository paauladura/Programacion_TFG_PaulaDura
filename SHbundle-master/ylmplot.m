function hh = ylmplot(l, m, shape,derivative)

% YLMPLOT creates a surface spherical harmonic for degree L and order M.
% The constituting Legendre function and cosine or sine of m*lambda
% are drawn as well.
%
% IN:
%    l ........ degree
%    m ........ order.  m>0 -> cos(m*lambda),  m<0 -> sin(m*lambda).
%    shape..... surface of projection
%               * 'plane':  2D plot with longitude/latitude grid
%               * 'sphere': projection on the unit sphere
%               * 'geoid':  up to 10% radial deviations of the sphere like 
%                           'Potsdamer Geoid'  
%               * 'orbital': Ylm as radial deviation of the origin
%                           like 'atom orbitals' in chemistry
%    derivative 
%               * 'dylm':  1. derivative of Ylm w.r.t. theta
%               * 'ddylm': 2. derivative of Ylm w.r.t. theta
%               * 'none'   (default) visualization of Ylm
% OUT:
%    h ........ vector of handles to 1) Ylm-axes, 2) Plm-axes, and 3) cos/sin-axes.
%
% USES:
%    PLM, uberall/ISINT
 
% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAGP, TU Munich
%    Matthias WEIGELT (MW),  GI, Uni Stuttgart 
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2020-12-01: MA, visualization of Ylm-derivatives
%    2014-12-02: MA, visualization on the sphere
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-29: MA, comments
%    2013-01-23: MA, input of plm in radian
%    2013-01-18: MA, internal colormap, xlabel, ylabel 
%    1995-07-06: NS, initial version
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

if ~isscalar(l) || ~isint(l), error(' L must be integer scalar'), end
if nargin == 1,   m = 0; end
if ~isscalar(m) || ~isint(m), error(' M must be integer scalar'), end
if abs(m) > l, error('m should be within [-l l]'), end
if nargin<3
    shape = 'plane';
end
if nargin<4
  derivative = 'none'
end
% create theta and lambda grid, depending on # zero crossings.
nlam   = min(400,(abs(m)+1)*100);
ntet   = min(200,(l-abs(m)+2)*50);
lambda = linspace(-pi,pi,nlam);			% [rad]
theta  = linspace(0,pi,ntet)';                  % [rad]
thetaDEG = theta*180/pi;
% create the functions.
[p,dp,ddp] = plm(l,abs(m),theta);			% Legendre function (standing) 
switch derivative
case 'ddylm'
  p = ddp;
  titel_plm = ['d^2P_{' num2str(l) ',' num2str(m) '}(cos \vartheta)/d\vartheta^2'];
  txt = '2x differentiated Surface Spherical Harmonic:';
case 'dylm'
  p = dp;
  titel_plm = ['dP_{' num2str(l) ',' num2str(m) '}(cos \vartheta)/d\vartheta'];
  txt = 'Differentiated Surface Spherical Harmonic:'
otherwise
  titel_plm = ['P_{' num2str(l) ',' num2str(m) '}(cos \vartheta)'];
  txt = 'Surface Spherical Harmonic:';
end
if m >= 0, c = cos(m*lambda); else c = sin(-m*lambda); end	% (lying)
y = p*c;					% surface spherical harmonic
lambdaDEG = lambda*180/pi;				% [deg]

% Initialize figure window
clf
set(gcf,'units','normalized')
n = length(colormap);
colormap(redblue1(n+rem(n,2)-1))		% length(colormap) is odd

% create the Ylm-axes
h(1) = axes('Position',[.3 .1  .6  .6 ]);
lamTitle = false;
switch shape 
case 'sphere'
    [Lam,The] = meshgrid(lambda,theta);
    [X,Y,Z] = sph2cart(Lam,pi/2-The,1);
    hs   = surf(X,Y,Z,y); shading interp;
    axis equal; axis off;
    lamTitle = true;
case 'geoid'
    MAX = 10*max(max(abs(y)));
    [Lam,The] = meshgrid(lambda,theta);
    [X,Y,Z] = sph2cart(Lam,pi/2-The,1+y/MAX);
    hs   = surf(X,Y,Z,y); shading interp;
    axis equal; axis off;
    lamTitle = true;
case 'orbital'
    [Lam,The] = meshgrid(lambda,theta);
    [X,Y,Z] = sph2cart(Lam,pi/2-The,y);
    hs   = surf(X,Y,Z,y); shading interp;
    axis equal; axis off;
    lamTitle = true;
otherwise
    %hs   = pcolor(lambda,thetaDEG,y);
    hs   = imagesc(lambdaDEG,thetaDEG,y);
    set(h(1),'YDir','reverse','XLim',[-180 180],'YLim',[0 180])
    set(h(1),'XTick',-180:60:180,'YTick',0:30:180,'Layer','top')
    
    shading flat
    grid
    xlabel('\lambda')
end

% create the Plm-axes
h(2) = axes('Position',[.1 .1  .15 .6 ]);
line(p,thetaDEG,'Color','r','LineWidth',3)
set(h(2),'YDir','reverse','YLim',[0 180],'YTick',0:30:180)
ylabel('\vartheta')
title(titel_plm)
grid

% create the cos/sin-axes
h(3) = axes('Position',[.3 .75 .6  .15]);
line(lambdaDEG,c,'Color','b','LineWidth',3)
set(h(3),'XLim',[-180 180],'XTickLabel',[],'XTick',-180:60:180)
if m >= 0, title(['cos(' num2str(m) '\lambda)']), else title(['sin(' num2str(m) '\lambda)']), end
if lamTitle == true,  xlabel('\lambda'); end
grid

if nargout == 1, hh = h; end
axes(h(1))					% current axes now
if exist('suptitle','file') == 2
   txt = [txt 'l = ' int2str(l)];
   txt = [txt ', m = ' int2str(m)];
   suptitle(txt)
end

end


function rgb = redblue1(n)

% redblue1(n) creates a colormap, ranging from dark blue via white to dark red.
%         
% Nico Sneeuw
% Munich, 31/08/94

% previously called seismic

if nargin == 0, n = size(get(gcf,'colormap'),1); end

m    = ceil(n/3);
top  = ones(m,1);
bot  = zeros(m,1);
up   = (0:m-1)'/m;
down = flipud(up);

r    = [bot; up; 1; top; down];
g    = [bot; up; 1; down; bot];
b    = [up; top; 1; down; bot];
rgb  = [r g b];

% rgb-map has size 4m+1 now. The central part will be extracted.

xlarge = 4*m+1-n;
xblue  = round(xlarge/2);
xred   = xlarge - xblue;
rgb([1:xblue 4*m-xred+2:4*m+1],:) = [];
end
