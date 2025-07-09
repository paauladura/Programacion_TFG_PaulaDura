function gamma = normg(phi, h, ellipsoid)

% function to calculate the normal gravity in latitude phi and heigh h
%
% IN:
%    phi ....... latitude in [rad]
%    h ......... heigh in [m]
%                ellipsoid: GRS80 or GRS67
%
% OUT:
%    gamma ..... the normal gravity in latitude phi and heigh h
%                also possible:
%                gamma=normg(phi,h)     GRS80 is default
%                gamma=normg(phi)       GRS80 is default and h=0

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), DoGE, UofC
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2000-01-31: NS, initial version
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

if nargin==2
   ellipsoid='GRS80';
end
if nargin==1
   ellipsoid='GRS80';
   h=0;
end

phi=phi(:);
if max(abs(phi))>pi/2
    warning('Is the latitude ''phi'' given in radian?')
end

if strcmp(ellipsoid,'GRS80')==1
   gamma_e=9.7803267715;
   gamma_p=9.8321863685;
   k=0.001931851353;
   e2=0.00669438002290;
   a=6378137.0;
   f=1/298.257222101;
   w=7.292115e-5;
   m=0.00344978600308;
   GM=3.986005e14;
elseif strcmp(ellipsoid,'GRS67')==1
   gamma_e=9.7803184558;
   gamma_p=9.8321772792;
   a=6378160.0;
   f=1/298.257222101;
   w=7.2921151467e-5;
   m=0.003358544730003;
   GM=3.98603e14;
   k=0.001931663383207;
   e2=0.006694605328561;
else
   error('Unkown ellipsoid')
end


% Berechnung von gamma in H=0;
gamma = gamma_e.*(1+k.*(sin(phi).^2))./(sqrt(1-e2.*(sin(phi).^2)));  %exakte Berechnung

%Berechnung von gamma in der Höhe H.
gamma=gamma.*(1-2/a.*(1+f+m-2*f*(sin(phi).^2)).*h + 3/a^2.*h.^2);
