function h = sctriplot(scmat,lmax,units)

% SCTRIPLOT plots the triangles of SH coefficients stored in SC-format. If
% stored in CS format it is converted to SC-format before plotting.
%
% sctriplot(scmat,lmax)
% h = sctriplot(scmat,lmax,units)
%
% INPUT
% scmat - Matrix of real SH coefficients in SC, CS or [l m Clm Slm] formats.
% lmax  - Maximum degree of spherical harmonic expansion [integer].
% units - Units to be labelled for the colorbar [string].
%
% OUTPUT
% Generates an image of the SC-formatted SH coefficients.
% h         -   A structure containing the handles for the following:
% h.fig     -   Handle of the figure.
% h.img     -   Handle of the grid image
% h.axis    -   Axis handle
% h.cbar    -   Colorbar handle
%-------------------------------------------------------------------------------

% 3 September 2007, Stuttgart
% Balaji Devaraju
%-------------------------------------------------------------------------------

if nargin == 2
    units = '';
elseif nargin == 1
    error('Insufficient input arguments')
end

[r,c] = size(scmat);

if (min(r,c) == lmax+1) && (r==c)
    scmat = cs2sc(scmat,0);
elseif (r==sum(1:lmax+1)) && (c==4)
    scmat = clm2sc(scmat);
elseif (r > c) && (r == (2*lmax + 1)) && (c == lmax+1)
    scmat = scmat';
elseif (min(r,c) ~= lmax+1) || (max(r,c) ~= (2*lmax + 1))
    error('Matrix neither confirms to SC-format, nor CS-format')
end

h.fig   = figure;
h.img   = imagesc(-lmax:lmax,0:lmax,log10(abs(scmat)));
h.axis  = gca();
h.cbar  = colorbar('Location', 'EastOutside');
pbaspect([2 1 1]);
set(get(h.cbar,'XLabel'),'String',units)          % changed
set(get(h.cbar,'YLabel'),'String','log10')        % changed
