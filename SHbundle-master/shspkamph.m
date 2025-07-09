function [amp, ph] = shspkamph(Klm,lmax,rc)

% SHSPKAMPH computes the amplitude and phase quantities of a given spherical
% harmonic spectrum.
%
% Klm 	= Alm * exp(1i * philm),
% Alm 	= abs(Klm);
% philm = arg(philm);
% where Klm is the complex spherical harmonic spectrum
% Alm is the amplitude part of the spectrum
% philm is the phase part of the spectrum
%
% IN: 
%    Klm ...... Spherical harmonic oefficients in C\S or S|C or [l m C S] formats
%    lmax ..... Maximum degree of spherical harmonic spectrum
%    rc ....... Flag to indicate whether the input Klm is in real or complex formats
% 				complex - 0		real - 1
%
% OUT:
%     amp ..... Alm in C\S format
%     ph ...... philm in C\S format [radians]
%
% REMARKS: 
%    All three inputs must be given!!!
%
% USES: 
%    real2cpxsh, cssc2clm, clm2sc, sc2cs 

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2021-04-12: MA, remove revision statement on deprecated function  
%    2013-04-08: BD, initial version
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

if rc == 1
	Klm = real2cpxsh(Klm,lmax);
elseif rc == 0
	[r, c] = size(Klm);
	if c == 4
		rklm = cssc2clm([Klm(:,1:2) real(Klm(:,3:4))],lmax);
		cklm = cssc2clm([Klm(:,1:2) imag(Klm(:,3:4))],lmax);
		rklm = sc2cs(clm2sc(rklm));
		cklm = sc2cs(clm2sc(cklm));
		Klm  = rklm + 1i*cklm;
	elseif c == (2*lmax  +  1) && r == (lmax+1)
		rklm = sc2cs(real(Klm));
		cklm = sc2cs(imag(Klm));
		Klm  = rklm + 1i*cklm;
	elseif c ~= r
		error('Input format of Klm is unknown')
	end
end

amp = abs(Klm);
ph 	= angle(Klm);

