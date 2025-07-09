%% EXAMPLE for the use of GSHSAG/GSHS_GRID - 2015-05-22
clear all;
clc;

if exist('gshs_grid','file') ~= 2,error('Please add path of ''SHBundle'''),   end; 
if exist('checkcoor','file') ~= 2,    error('Please add path of ''uberall'''),    end; 

constants; % load constants from UBERALL/constants

% parse ICGEM-file up to max_lm = 10
[gsm, lmax, lmin, info] = parse_icgem('./data/example.icgem', 'max_lm', 10); 
% convert l,m-indexed list into /S|C\-format
[field_sc, lmax]        = clm2sc(gsm, 'max_lm', lmax); 

% prepare a longitude/latitude grid
lam = -179:1:180; lRAD = lam * pi / 180;
phi = -90:1:90;   pRAD = phi * pi / 180;
[Lam, Phi] = meshgrid(lam, phi);


quantity = {'geoid','dg','tr','trr','slope','potential'}
%% calculation
figure
for ii = 1:6
    tic; 
    F   = gshs_grid(field_sc, lRAD, pRAD, ae, 'GM', GM, 'height', 0, 'max_lm', lmax, 'quant', quantity{ii}, 'sub_wgs84', true, 'curvature', false, 'waitbar', true, 'legendre', 'plm'); 
    toc;
    subplot(2, 3, ii);
    imagesc(lam, phi, F); shading flat; 
    title(['GSHS\_GRID: ' quantity{ii} ]); 
    axis equal; axis tight; colorbar
end

tic
Fc  = gshs_grid(field_sc, lRAD, pRAD, ae, 'GM', GM, 'height', 0, 'max_lm', lmax, 'quant', 'potential', 'sub_wgs84', true, 'curvature', true, 'waitbar', true, 'legendre', 'plm');
toc

figure;
subplot(2, 1, 1);
imagesc(lam, phi, F); shading flat; 
xlabel('longitude'); ylabel('latitude'); title('GSHS\_GRID: GRS80 reduced potential, [m^2/s^2]'); 
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[-90:30:90]); 
axis equal; colorbar



subplot(2, 1, 2);
imagesc(lam, phi, log10(abs(Fc))); shading flat; 
xlabel('longitude'); ylabel('latitude'); title('log10 of  curvature of GRS80 reduced potential'); 
set(gca,'xTick',[-180:30:180]); set(gca,'yTick',[-90:30:90]); 
axis equal; colorbar