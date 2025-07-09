%% EXAMPLE for the use of GSHS_                       - 2015-05-22
%  removed the demonstration of depricated functions  - 2021-04-12
clear all;
close all;
clc;

addpath('..');            % path to 'SHbundle'
addpath('../../uberall'); % path to 'uberall'

[gsm, lmax, lmin, info] = parse_icgem('./data/example.icgem', 'max_lm', 10); % read model
[field, lmax] = clm2sc(gsm, 'max_lm', lmax); % convert format of coeffs to /S|C\

%% computation
[V_pot, theRAD, lamRAD] = gshs_(field, 'quant', 'potential', 'grid', 'neumann', 'gridsize', lmax, 'height', 0, 'sub_wgs84', true);
[V_dg] = gshs_(field, 'quant', 'dg', 'grid', 'neumann', 'gridsize', lmax, 'height', 0, 'sub_wgs84', 1);
[V_rr] = gshs_(field, 'quant', 'trr', 'grid', 'neumann', 'gridsize', lmax, 'height', 0, 'sub_wgs84', 1);
% % also possible (using default values):
% [V_pot, theRAD, lamRAD] = gshs_(field, 'grid', 'neumann');
% [V_dg] = gshs_(field, 'quantity', 'dg', 'grid', 'neumann');
% [V_rr] = gshs_(field, 'quantity', 'trr', 'grid', 'neumann');

%% output
figure;
subplot(2, 2, 1);
imagesc(lamRAD*180/pi, theRAD*180/pi, V_pot); title('potential'); c = colorbar; ylabel(c, '[m^2/s^2]');
subplot(2, 2, 2);
imagesc(lamRAD*180/pi, theRAD*180/pi, V_dg); title('gravity anomaly'); c = colorbar; ylabel(c, '[m]');
subplot(2, 2, 3);
imagesc(lamRAD*180/pi, theRAD*180/pi, V_rr); title('2nd radial derivative'); c = colorbar; ylabel(c, '[1/s^2]');
