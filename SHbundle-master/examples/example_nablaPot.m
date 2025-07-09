%% EXAMPLE for the use of NABLAPOT - 2014-10-27
clear all;
close all;
clc;

addpath('..');            % path to 'SHbundle'
addpath('../../uberall'); % path to 'uberall'

[gsm, lmax, lmin, info] = parse_icgem('./data/example.icgem', 'max_lm', 10); % read ICGEM-model
[field, lmax] = clm2sc(gsm, 'max_lm', lmax);                                 % convert from clm to /S|C\ format

% (extremly drifting) circular orbit:
x = (0:.01:10*pi)'; inc = -85;  omasc = x/15;
r = 7200000; z = 0*x;
X = multmatvek([cosd(z+inc), z, -sind(z+inc),z z+1, z, sind(z+inc), z, cosd(z+inc)], [r.*cos(x), r.*sin(x), z]);
X = multmatvek([cosd(omasc), sind(omasc),z, -sind(omasc),  cosd(omasc), z, z, z, z+1], X);
[lamRAD, phiRAD, r] = cart2sph(X(:,1), X(:,2), X(:,3));

% computation 
tic; [Tr, Tth, Tlam, T] = nablaPot(field, lamRAD, pi/2 - phiRAD, r); toc;

% output
figure;
subplot(311); plot(Tr); title('radial derivative T_r'); 
subplot(312); plot(Tth); title('derivative in co-latitude direction: T_\theta'); 
subplot(313); plot(Tlam); title('derivative in longitude direction: T_\lambda'); 