function [posA, posB] = prpcoord(crdfname)

% PRPCOORD prepares the cartesian co-ordinates of the GRACE satellite
% orbits for the calculation of the normal matrices and error matrix
% simulation.
%
% [posA,posB] = prpcoord(crdfname)
% INPUT
% crdfname  -   Filename of the MAT file where the co-ordinates of each
%               month are stored.
% 
% OUTPUT
% posA,posB -   Positions of GRACE A&B satellite with lam(0-360),
%               theta(0-180), radius(m).
%--------------------------------------------------------------------------

% Created on: 29 August 2007, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

mats = load(crdfname);
field = fieldnames(mats);

[lam, th, r] = cart2sph(mats.(field{1})(:,3), mats.(field{1})(:,4), mats.(field{1})(:,5));
cnt = find(lam<0);
lam(cnt) = lam(cnt) + 2*pi;
th = pi/2 - th;
posA = [lam,th,r];

[lam, th, r] = cart2sph(mats.(field{2})(:,3), mats.(field{2})(:,4), mats.(field{2})(:,5));
cnt = find(lam<0);
lam(cnt) = lam(cnt) + 2*pi;
th = pi/2 - th;
posB = [lam,th,r];